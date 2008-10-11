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
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FeatureRealizer;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.impl.FeatureImpl;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.DummySymbolList;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;

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
 * @since 1.2
 * @author Thomas Down
 */

class RawAnnotationSequence
  extends
    Unchangeable
  implements
    DASSequenceI
{
    private URL dataSourceURL;
    private String seqID;
    private FeatureRealizer featureRealizer = FeatureImpl.DEFAULT;
    private SymbolList nullSymbols;
    private DASSequenceDB dummyDB;
    private DASFeatureSet features;

    {
        nullSymbols = new DummySymbolList(DNATools.getDNA(), 2000000000);
    }


    RawAnnotationSequence(DASSequenceDB dummyDB, URL dataSourceURL, String seqID)
        throws BioException, IllegalIDException
    {
        this.dummyDB = dummyDB;
        this.dataSourceURL = dataSourceURL;
        this.seqID = seqID;
        this.features = new DASFeatureSet(this, dataSourceURL, seqID);
    }

    URL getDataSourceURL() {
        return dataSourceURL;
    }

    public DASSequenceDB getParentDB() {
        return dummyDB;
    }


    private int registerLocalFeatureFetchers(Object regKey) {
        features.registerFeatureFetcher(regKey);

        return 1;
    }

    private int registerLocalFeatureFetchers(Location l, Object regKey) {
        features.registerFeatureFetcher(l, regKey);

        return 1;
    }

    int registerFeatureFetchers(Object regKey) throws BioException {
        return registerLocalFeatureFetchers(regKey);
    }

    int registerFeatureFetchers(Location l, Object regKey) throws BioException {
        return registerLocalFeatureFetchers(l, regKey);
    }

    //
    // SymbolList stuff
    //

    public Alphabet getAlphabet() {
        return nullSymbols.getAlphabet();
    }

    public Iterator iterator() {
        return nullSymbols.iterator();
    }

    public int length() {
        return nullSymbols.length();
    }

    public String seqString() {
        return nullSymbols.seqString();
    }

    public String subStr(int start, int end) {
        return nullSymbols.subStr(start, end);
    }

    public SymbolList subList(int start, int end) {
        return nullSymbols.subList(start, end);
    }

    public Symbol symbolAt(int pos) {
        return nullSymbols.symbolAt(pos);
    }

    public List toList() {
        return nullSymbols.toList();
    }

    public void edit(Edit e)
        throws ChangeVetoException
    {
        throw new ChangeVetoException("/You/ try implementing read-write DAS");
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
          throw new BioRuntimeException("Couldn't create features iterator",be);
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

            FeatureFilter structureMembershipFilter = new FeatureFilter.ByClass(ComponentFeature.class);

            if (FilterUtils.areProperSubset(ff, structureMembershipFilter)) {
                return FeatureHolder.EMPTY_FEATURE_HOLDER;
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
        return featureRealizer.realizeFeature(this, dest, temp);
    }

    //
    // Annotatable stuff
    //

    public Annotation getAnnotation() {
        return Annotation.EMPTY_ANNOTATION;
    }
}
