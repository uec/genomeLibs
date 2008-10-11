/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

package org.biojava.bio.program.das;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.SAXParserFactory;

import org.biojava.bio.Annotatable;
import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FeatureRealizer;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.LazyFeatureHolder;
import org.biojava.bio.seq.MergeFeatureHolder;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.impl.AssembledSymbolList;
import org.biojava.bio.seq.impl.FeatureImpl;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.SeqIOAdapter;
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.AbstractChangeable;
import org.biojava.utils.ChangeEvent;
import org.biojava.utils.ChangeForwarder;
import org.biojava.utils.ChangeSupport;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.cache.CacheReference;
import org.xml.sax.XMLReader;

/**
 * Sequence reflecting a DAS reference sequence, possibly
 * decorated with one of more annotation sets.
 *
 * <p>
 * This is an first-pass implementation.  In future, I hope
 * to add query optimization for better performance on large
 * sequences, and pluggable transducers to parameterize the
 * creation of BioJava features.
 * </p>
 *
 * @since 1.1
 * @author Thomas Down
 * @author Matthew Pocock
 * @author David Huen
 */

public class DASSequence
  extends
    AbstractChangeable
  implements
    DASSequenceI
{
    /**
     * Change type which indicates that the set of annotation servers used
     * by this DASSequence has been changed. This extends Feature.FEATURES as
     * the addition and removal of annotation servers adds and removes features.
     */

    public static final ChangeType ANNOTATIONS = new ChangeType(
            "Annotation sets have been added or removed from the DAS sequence",
            "org.biojava.bio.program.das.DASSequence",
            "ANNOTATIONS",
            Feature.FEATURES
    );

    public static final String PROPERTY_ANNOTATIONSERVER = "org.biojava.bio.program.das.annotation_server";
    public static final String PROPERTY_FEATUREID = "org.biojava.bio.program.das.feature_id";
    public static final String PROPERTY_FEATURELABEL = "org.biojava.bio.program.das.feature_label";
    public static final String PROPERTY_LINKS = "org.biojava.bio.program.das.links";
    public static final String PROPERTY_SEQUENCEVERSION = "org.biojava.bio.program.das.sequence_version";

    public static final int SIZE_THRESHOLD = 500000;

    private static final int SYMBOL_TILE_THRESHOLD = 100000;
    private static final int SYMBOL_TILE_SIZE = 20000;

    protected transient ChangeForwarder annotationForwarder;

    private DASSequenceDB parentdb;
    private URL dataSourceURL;
    private String seqID;
    private String version = null;
    private FeatureRealizer featureRealizer = FeatureImpl.DEFAULT;
    private FeatureRequestManager.Ticket structureTicket;

    private CacheReference refSymbols;
    private int length = -1;

    private Map featureSets;
    private FeatureHolder structure;
    private MergeFeatureHolder features;

    {
        featureSets = new HashMap();
        features = new MergeFeatureHolder();
    }

    private class StructureWrapper extends LazyFeatureHolder {
        public StructureWrapper() {
            super(new FeatureFilter.And(
                    FeatureFilter.top_level,
                    new FeatureFilter.ByClass(ComponentFeature.class)
            ));
        }

        protected FeatureHolder createFeatureHolder() {
            try {
                return getStructure();
            } catch (BioException ex) {
                throw new BioRuntimeException("Error fetching structure",ex);
            }
        }
    }

    DASSequence(DASSequenceDB db, URL dataSourceURL, String seqID, Set dataSources)
        throws BioException, IllegalIDException
    {
        this.parentdb = db;
        this.dataSourceURL = dataSourceURL;
        this.seqID = seqID;

        //
        // Check for deep structure.  This also checks that the sequence
        // really exists, and hopefully picks up the length along the way.
        //

        SeqIOListener listener = new SkeletonListener();
        FeatureRequestManager frm = getParentDB().getFeatureRequestManager();
        this.structureTicket = frm.requestFeatures(getDataSourceURL(), seqID, listener, null, "component");
    try {
        features.addFeatureHolder(new StructureWrapper());
    } catch (ChangeVetoException cve) {
                throw new BioError(cve);
        }

        //
        // Pick up some annotations
        //

    for (Iterator dsi = dataSources.iterator(); dsi.hasNext(); ) {
        URL annoURL = (URL) dsi.next();

        FeatureHolder newFeatureSet = new DASFeatureSet(this, annoURL, seqID);
        featureSets.put(annoURL, newFeatureSet);
        try {
            features.addFeatureHolder(newFeatureSet);
        } catch (ChangeVetoException cve) {
            throw new BioError(cve);
        }
    }
    }

    private class SkeletonListener extends SeqIOAdapter {
        private SimpleFeatureHolder structureF;

        public void startSequence() {
            structureF = new SimpleFeatureHolder();
        }

        public void endSequence() {
            structure = structureF;
            // if (structure.countFeatures() > 0) {
                // try {
                    // features.addFeatureHolder(structure, new FeatureFilter.ByClass(ComponentFeature.class));
                // } catch (ChangeVetoException cve) {
                //    throw new BioError(cve);
                // }
            // }
        }

        public void addSequenceProperty(Object key, Object value)
            throws ParseException
        {
            try {
                if (key.equals("sequence.start")) {
                    int start = Integer.parseInt(value.toString());
                    if (start != 1) {
                        throw new ParseException("Server doesn't think sequence starts at 1.  Wierd.");
                    }
                } else if (key.equals("sequence.stop")) {
                    length = Integer.parseInt(value.toString());
                } else if (key.equals("sequence.version")) {
                    version = value.toString();
                }
            } catch (NumberFormatException ex) {
                throw new ParseException(ex, "Expect numbers for segment start and stop");
            }
        }

        public void startFeature(Feature.Template temp)
            throws ParseException
        {
            if (temp instanceof ComponentFeature.Template) {

                try {
                    ComponentFeature.Template ctemp = (ComponentFeature.Template) temp;
                    ComponentFeature cf = new DASComponentFeature(DASSequence.this,
                                                                  ctemp);
                    structureF.addFeature(cf);

                    length = Math.max(length, ctemp.location.getMax());
                } catch (BioException ex) {
                    throw new ParseException(ex, "Error instantiating DASComponent");
                } catch (ChangeVetoException ex) {
                    throw new BioError("Immutable FeatureHolder when trying to build structure",ex);
                }
            } else {
                // Server seems not to honour category=
                // This hurts performance, but we can just elide the unwanted
                // features on the client side.
            }
        }
    }

    URL getDataSourceURL() {
        return dataSourceURL;
    }

    public DASSequenceDB getParentDB() {
        return parentdb;
    }

    FeatureHolder getStructure() throws BioException {
        if(!this.structureTicket.isFetched()) {
            this.structureTicket.doFetch();
        }

        return this.structure;
    }

    private void _addAnnotationSource(URL dataSourceURL)
        throws BioException, ChangeVetoException
    {
        FeatureHolder structure = getStructure();
        for (Iterator i = structure.features(); i.hasNext(); ) {
            DASComponentFeature dcf = (DASComponentFeature) i.next();
            DASSequence seq = dcf.getSequenceLazy();
            if (seq != null) {
                seq.addAnnotationSource(dataSourceURL);
            }
        }

        FeatureHolder fs = new DASFeatureSet(this, dataSourceURL, this.seqID);
        featureSets.put(dataSourceURL, fs);
        features.addFeatureHolder(fs);
    }

    public Set dataSourceURLs() {
      return Collections.unmodifiableSet(featureSets.keySet());
    }

    public void addAnnotationSource(URL dataSourceURL)
        throws BioException, ChangeVetoException
    {
      if(!featureSets.containsKey(dataSourceURL)) {
        if (!hasListeners()) {
          _addAnnotationSource(dataSourceURL);
        } else {
          ChangeSupport changeSupport = getChangeSupport(ANNOTATIONS);
          synchronized (changeSupport) {
            ChangeEvent ce = new ChangeEvent(
            this,
            ANNOTATIONS,
            null,
            null
            ) ;
            changeSupport.firePreChangeEvent(ce);
            _addAnnotationSource(dataSourceURL);
            changeSupport.firePostChangeEvent(ce);
          }
        }
      }
    }

    protected ChangeSupport getChangeSupport(ChangeType ct){
      ChangeSupport cs = super.getChangeSupport(ct);

      if(annotationForwarder == null &&
        (ct == null || ct == Annotatable.ANNOTATION)){
        annotationForwarder = new ChangeForwarder.Retyper(this, cs, Annotation.PROPERTY);
        getAnnotation().addChangeListener(
            annotationForwarder,
            Annotatable.ANNOTATION);
      }
      return cs;
    }

    private void _removeAnnotationSource(URL dataSourceURL)
        throws ChangeVetoException, BioException
    {
        FeatureHolder structure = getStructure();
        FeatureHolder fh = (FeatureHolder) featureSets.get(dataSourceURL);
        if (fh != null) {
            for (Iterator i = structure.features(); i.hasNext(); ) {
                DASComponentFeature dcf = (DASComponentFeature) i.next();
                DASSequence seq = dcf.getSequenceLazy();
                if (seq != null) {
                    seq.removeAnnotationSource(dataSourceURL);
                }
            }

            features.removeFeatureHolder(fh);
            featureSets.remove(dataSourceURL);
        }
    }

    public void removeAnnotationSource(URL dataSourceURL)
        throws ChangeVetoException, BioException
    {
      if (featureSets.containsKey(dataSourceURL)) {
        if (!hasListeners()) {
          _removeAnnotationSource(dataSourceURL);
        } else {
          ChangeSupport changeSupport = getChangeSupport(ANNOTATIONS);
          synchronized (changeSupport) {
                      ChangeEvent ce = new ChangeEvent(
                                                     this,
                                                     ANNOTATIONS,
                                                     null,
                                                     null
                                                     ) ;
                      changeSupport.firePreChangeEvent(ce);
          _removeAnnotationSource(dataSourceURL);
          changeSupport.firePostChangeEvent(ce);
          }
        }
      }
    }

    private int registerLocalFeatureFetchers(Object regKey) {
        // System.err.println(getName() + ": registerLocalFeatureFetchers()");

        for (Iterator i = featureSets.values().iterator(); i.hasNext(); ) {
            DASFeatureSet dfs = (DASFeatureSet) i.next();
            dfs.registerFeatureFetcher(regKey);
        }

        return featureSets.size();
    }

    private int registerLocalFeatureFetchers(Location l, Object regKey) {
        // System.err.println(getName() + ": registerLocalFeatureFetchers(" + l.toString() + ")");

        for (Iterator i = featureSets.values().iterator(); i.hasNext(); ) {
            DASFeatureSet dfs = (DASFeatureSet) i.next();
            dfs.registerFeatureFetcher(l, regKey);
        }

        return featureSets.size();
    }

    int registerFeatureFetchers(Object regKey) throws BioException {
        // System.err.println(getName() + ": registerFeatureFetchers()");

        for (Iterator i = featureSets.values().iterator(); i.hasNext(); ) {
            DASFeatureSet dfs = (DASFeatureSet) i.next();
            dfs.registerFeatureFetcher(regKey);
        }

        int num = featureSets.size();

        FeatureHolder structure = getStructure();
        if (length() < SIZE_THRESHOLD && structure.countFeatures() > 0) {
            List sequences = new ArrayList();
            for (Iterator fi = structure.features(); fi.hasNext(); ) {
                ComponentFeature cf = (ComponentFeature) fi.next();
                DASSequence cseq = (DASSequence) cf.getComponentSequence();
                sequences.add(cseq);
            }

            for (Iterator si = sequences.iterator(); si.hasNext(); ) {
                DASSequence cseq = (DASSequence) si.next();
                num += cseq.registerFeatureFetchers(regKey);
            }
        }

        return num;
    }

    int registerFeatureFetchers(Location l, Object regKey) throws BioException {
        // System.err.println(getName() + ": registerFeatureFetchers(" + l.toString() + ")");

        for (Iterator i = featureSets.values().iterator(); i.hasNext(); ) {
            DASFeatureSet dfs = (DASFeatureSet) i.next();
            dfs.registerFeatureFetcher(l, regKey);
        }

        int num = featureSets.size();

        FeatureHolder structure = getStructure();
        if (structure.countFeatures() > 0) {
            FeatureHolder componentsBelow = structure.filter(new FeatureFilter.OverlapsLocation(l), false);

            Map sequencesToRegions = new HashMap();

            for (Iterator fi = componentsBelow.features(); fi.hasNext(); ) {
                ComponentFeature cf = (ComponentFeature) fi.next();
                DASSequence cseq = (DASSequence) cf.getComponentSequence();
                if (l.contains(cf.getLocation())) {
                    sequencesToRegions.put(cseq, null);
                } else {
                    Location partNeeded = l.intersection(cf.getLocation());
                    if (cf.getStrand() == StrandedFeature.POSITIVE) {
                        partNeeded = partNeeded.translate(cf.getComponentLocation().getMin() - cf.getLocation().getMin());
                        sequencesToRegions.put(cseq, partNeeded);
                    } else {
                        sequencesToRegions.put(cseq, null);
                    }
                }
            }

            for (Iterator sri = sequencesToRegions.entrySet().iterator(); sri.hasNext(); ) {
                Map.Entry srme = (Map.Entry) sri.next();
                DASSequence cseq = (DASSequence) srme.getKey();
                Location partNeeded = (Location) srme.getValue();
                if (partNeeded != null) {
                    num += cseq.registerFeatureFetchers(partNeeded, regKey);
                } else {
                    num += cseq.registerFeatureFetchers(regKey);
                }
            }
        }

        return num;
    }

    //
    // SymbolList stuff
    //

    public Alphabet getAlphabet() {
        try {
            return getSymbols().getAlphabet();
        } catch (BioException be) {
            throw new BioRuntimeException(be);
        }
    }

    public Iterator iterator() {
      try {
        return getSymbols().iterator();
      } catch (BioException be) {
        throw new BioRuntimeException("Can't iterate over symbols", be);
      }
    }

    public int length() {
        try {
            if (length < 0 && !structureTicket.isFetched()) {
                // Hope that the length is set when we get the structure

                structureTicket.doFetch();
            }

            if (length < 0) {
                throw new BioError("Assertion Failure: structure fetch didn't get length");
            }

            return length;
        } catch (BioException be) {
            throw new BioRuntimeException("Can't calculate length",be);
        }
    }

    public String seqString() {
        try {
            return getSymbols().seqString();
        } catch (BioException be) {
            throw new BioRuntimeException("Can't create seqString",be);
        }
    }

    public String subStr(int start, int end) {
        try {
            return getSymbols().subStr(start, end);
        } catch (BioException be) {
            throw new BioRuntimeException("Can't create substring",be);
        }
    }

    public SymbolList subList(int start, int end) {
        try {
            return getSymbols().subList(start, end);
        } catch (BioException be) {
            throw new BioRuntimeException("Can't create subList",be);
        }
    }

    public Symbol symbolAt(int pos) {
        try {
            return getSymbols().symbolAt(pos);
        } catch (BioException be) {
            throw new BioRuntimeException("Can't fetch symbol",be);
        }
    }

    public List toList() {
        try {
            return getSymbols().toList();
        } catch (BioException be) {
            throw new BioRuntimeException("Can't create list",be);
        }
    }

    public void edit(Edit e)
        throws ChangeVetoException
    {
        throw new ChangeVetoException("/You/ try implementing read-write DAS");
    }

    //
    // DNA fetching stuff
    //

    protected SymbolList getSymbols() throws BioException {
        SymbolList sl = null;
        if (refSymbols != null) {
            sl = (SymbolList) refSymbols.get();
        }

        if (sl == null) {
            FeatureHolder structure = getStructure();
            if (structure.countFeatures() == 0) {
                if (length() > SYMBOL_TILE_THRESHOLD) {
                    AssembledSymbolList asl = new AssembledSymbolList();
                    int tileStart = 1;
                    while (tileStart <= length()) {
                        int tileEnd = Math.min(tileStart + SYMBOL_TILE_SIZE - 1, length());
                        if ((length() - tileEnd) < 1000) {
                            // snap the tiny bit onto the end of the last tile
                            tileEnd = length();
                        }
                        Location l = new RangeLocation(tileStart, tileEnd);
                        SymbolList symbols = new DASRawSymbolList(this,
                                                                  new Segment(getName(), tileStart, tileEnd));
                        asl.putComponent(l, symbols);

                        tileStart = tileEnd + 1;
                    }

                    sl = asl;
                } else {
                    sl = new DASRawSymbolList(this, new Segment(getName()));
                }
            } else {
                //
                // VERY, VERY, naive approach to identifying canonical components.  FIXME
                //

                Location coverage = Location.empty;
                AssembledSymbolList asl = new AssembledSymbolList();
                asl.setLength(length);

                for (Iterator i = structure.features(); i.hasNext(); ) {
                    ComponentFeature cf = (ComponentFeature) i.next();
                    Location loc = cf.getLocation();
                    if (LocationTools.overlaps(loc, coverage)) {
                        // Assume that this is non-cannonical.  Maybe not a great thing to
                        // do, but...
                    } else {
                        asl.putComponent(cf);
                        coverage = LocationTools.union(coverage, loc);
                    }
                }

                sl = asl;
            }

            refSymbols = parentdb.getSymbolsCache().makeReference(sl);
        }

        return sl;
    }

    //
    // Identification stuff
    //

    public String getName() {
        return seqID;
    }

    public String getURN() {
        try {
            return new URL(getDataSourceURL(), "?ref=" + seqID).toString();
        } catch (MalformedURLException ex) {
            throw new BioRuntimeException(ex);
        }
    }

    //
    // FeatureHolder stuff
    //

    public Iterator features() {
      try {
        registerFeatureFetchers(null);
        return features.features();
      } catch (BioException be) {
        throw new BioRuntimeException("Couldn't create features iterator", be);
      }
    }

    public boolean containsFeature(Feature f) {
        return features.containsFeature(f);
    }

    public FeatureHolder filter(FeatureFilter ff) {
        return filter(ff, !FilterUtils.areProperSubset(ff, FeatureFilter.top_level));
    }

    public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
        try {
            //
            // We optimise for the case of just wanting `structural' features,
            // which improves the scalability of the Dazzle server (and probably
            // other applications, too)
            //

            FeatureHolder structure = getStructure();
            FeatureFilter structureMembershipFilter = new FeatureFilter.ByClass(ComponentFeature.class);

            if (FilterUtils.areProperSubset(ff, structureMembershipFilter)) {
                if (recurse) {
                    for (Iterator fi = structure.features(); fi.hasNext(); ) {
                        //ComponentFeature cf = (ComponentFeature) fi.next();
                        //Sequence cseq = cf.getComponentSequence();
                        // Just ensuring that the sequence is instantiated should be sufficient.
                    }
                }
                return structure.filter(ff, recurse);
            }

            //
            // Otherwise they want /real/ features, I'm afraid...
            //

            Location ffl = FilterUtils.extractOverlappingLocation(ff);
            if (recurse) {
                int numComponents = 1;
                if (ffl != null) {
                    numComponents = registerFeatureFetchers(ffl, ff);
                } else {
                    numComponents = registerFeatureFetchers(ff);
                }
                getParentDB().ensureFeaturesCacheCapacity(numComponents * 3);
            } else {
                if (ffl != null) {
                    registerLocalFeatureFetchers(ffl, ff);
                } else {
                    // new Exception("Hmmm, this is dangerous...").printStackTrace();
                    registerLocalFeatureFetchers(ff);
                }
            }

            return features.filter(ff, recurse);
        } catch (BioException be) {
            throw new BioRuntimeException("Can't filter",be);
        }
    }

    public FeatureFilter getSchema() {
        return features.getSchema();
    }

    public int countFeatures() {
        return features.countFeatures();
    }

    public Feature createFeature(Feature.Template temp)
        throws ChangeVetoException
    {
        throw new ChangeVetoException("Can't create features on DAS sequences.");
    }

    public void removeFeature(Feature f)
        throws ChangeVetoException
    {
        throw new ChangeVetoException("Can't remove features from DAS sequences.");
    }

    //
    // Feature realization stuff
    //

    public Feature realizeFeature(FeatureHolder dest,
                                  Feature.Template temp)
        throws BioException
    {
        // workaround for features with bad locations - usualy starting at 0
        // e.g. enxembl830 has a cytogenetic band starting at 0
        if(temp.location.getMin() < 1 || temp.location.getMax() > length()) {
            temp.location = LocationTools.intersection(temp.location,
                                                       new RangeLocation(1, length()) );
        }

        return featureRealizer.realizeFeature(this, dest, temp);
    }

    //
    // Annotatable stuff
    //

    public Annotation getAnnotation() {
        try {
            Annotation anno = new SmallAnnotation();
            anno.setProperty(PROPERTY_SEQUENCEVERSION, version);
            return anno;
        } catch (ChangeVetoException ex) {
            throw new BioError("Expected to be able to modify annotation");
        }
    }

    //
    // Utility method to turn of the awkward bits of Xerces-J
    //

    static DocumentBuilder nonvalidatingParser() {
        try {
            DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
            dbf.setValidating(false);
            return dbf.newDocumentBuilder();
        } catch (Exception ex) {
            throw new BioError(ex);
        }
    }

    static XMLReader nonvalidatingSAXParser() {
        try {
            SAXParserFactory spf = SAXParserFactory.newInstance();
            spf.setValidating(false);
            spf.setNamespaceAware(true);
            return spf.newSAXParser().getXMLReader();
        } catch (Exception ex) {
            throw new BioError(ex);
        }
    }
}
