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

import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.bio.Annotation;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.SmallAnnotation;
import org.biojava.bio.program.xff.XFFFeatureSetHandler;
import org.biojava.bio.seq.ComponentFeature;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.MergeFeatureHolder;
import org.biojava.bio.seq.RealizingFeatureHolder;
import org.biojava.bio.seq.SimpleFeatureHolder;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.seq.io.SeqIOAdapter;
import org.biojava.bio.seq.io.SeqIOListener;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;
import org.biojava.utils.cache.CacheReference;

/**
 * FeatureHolder reflecting features provided by a DAS annotation
 * server.
 *
 * @since 1.1
 * @author Thomas Down
 * @author Matthew Pocock
 */

class DASFeatureSet
  extends
    Unchangeable
  implements
    FeatureHolder
{
    private FeatureRequestManager.Ticket[] featureTickets;
    private Location[]                     tiles;
    private CacheReference[]               tileFeatures;
    private SimpleFeatureHolder            unrulyFeatures = new SimpleFeatureHolder();
    private FeatureHolder                  allFeatures;

    private Map                            typesMap;
    private FeatureRequestManager.Ticket   typesTicket;
    private FeatureFilter                  allTypesFilter;

    private DASSequenceI                   refSequence;
    private URL                            dataSource;
    private String                         sourceID;

    DASFeatureSet(DASSequenceI seq, URL ds, String id)
        throws BioException
    {
        refSequence = seq;
        dataSource = ds;
        sourceID = id;
    }

    private final static int TILE_THRESHOLD_LENGTH   = 1000000;
    private final static int TILE_THRESHOLD_COUNT    = 2000;
    private final static int TILE_SIZE               = 100000;

    private Map getTypesMap()
    {
        if (typesMap == null) {
            if (typesTicket == null) {
                FeatureRequestManager frm = refSequence.getParentDB().getFeatureRequestManager();
                typesTicket = frm.requestTypes(dataSource,
                                               new Segment(refSequence.getName()),
                                               new DASTypesPopulator());
            }

            try {
                typesTicket.doFetch();
            } catch (Exception ex) {
                // Really evil hack-around for wormbase.

                ex.printStackTrace();
                // throw new BioRuntimeException(ex, "Error fetching types");
                try {
                    Set allTypes = DAS.getTypes(dataSource);
                    typesMap = new HashMap();
                    for (Iterator i = allTypes.iterator(); i.hasNext(); ) {
                        String type = (String) i.next();
                        typesMap.put(type, null);
                    }
                } catch (BioException ex2) {
                    throw new BioRuntimeException( "Types command isn't working AT ALL!", ex2);
                }
            }
        }
        if (typesMap == null) {
            throw new BioError("Assertion failure: types fetch hasn't happened yet");
        }

        return typesMap;
    }

    private FeatureFilter getAllTypesFilter() {
        if (allTypesFilter == null) {
            allTypesFilter = FeatureFilter.all;

            if (refSequence.length() > TILE_THRESHOLD_LENGTH) {
                Map typesMap = getTypesMap();
                for (Iterator ti = typesMap.keySet().iterator(); ti.hasNext(); ) {
                    String type = (String) ti.next();
                    FeatureFilter typeFilter = new FeatureFilter.ByType(type);
                    if (allTypesFilter == FeatureFilter.all) {
                        allTypesFilter = typeFilter;
                    } else {
                        allTypesFilter = new FeatureFilter.Or(allTypesFilter, typeFilter);
                    }
                }
            }
        }

        return allTypesFilter;
    }


    private Location[] getTiles() {
        if (tiles == null) {
            boolean doTiling = false;
            int seqLength = refSequence.length();

            if (seqLength > TILE_THRESHOLD_LENGTH) {
                // System.err.print("*** Considering tiling...");
                Map types = getTypesMap();

                int totalCount = 0;
                for (Iterator ti = types.values().iterator(); ti.hasNext(); ) {
                    Integer count = (Integer) ti.next();
                    if (count != null) {
                        totalCount += count.intValue();
                    } else {
                        doTiling = true;
                    }
                }

                if (doTiling) {
                    // System.err.println("yes (unknown total)");
                } else {
                    doTiling = (totalCount > TILE_THRESHOLD_COUNT);
                    //  if (doTiling) {
//  			System.err.println("yes (" + totalCount + ")");
//  		    } else {
//  			System.err.println("no.");
//  		    }
                }
            }

            if (doTiling)
            {
                int numTiles = (int) Math.ceil(1.0 * seqLength / TILE_SIZE);
                tiles = new Location[numTiles];
                featureTickets = new FeatureRequestManager.Ticket[numTiles];
                tileFeatures = new CacheReference[numTiles];
                for (int i = 0; i < numTiles; ++i) {
                    tiles[i] = new RangeLocation(i * TILE_SIZE + 1, Math.min((i + 1) * TILE_SIZE + 1, seqLength));
                }
            } else {
                tiles = new Location[1];
                tiles[0] = new RangeLocation(1, seqLength);
            }

            featureTickets = new FeatureRequestManager.Ticket[tiles.length];
            tileFeatures = new CacheReference[tiles.length];
        }

        return tiles;
    }

    private void registerFeatureFetcher(int tileNum, Object regKey) {
        Location[] tiles = getTiles();

        if (tileFeatures[tileNum] == null || tileFeatures[tileNum].get() == null) {
            if (featureTickets[tileNum] == null) {
                SeqIOListener listener = new DASFeatureSetPopulator(tileNum);
                FeatureRequestManager frm = refSequence.getParentDB().getFeatureRequestManager();
                if (tiles.length > 1) {
                    featureTickets[tileNum] = frm.requestFeatures(dataSource,
                                                                  sourceID,
                                                                  listener,
                                                                  tiles[tileNum]);
                } else {
                    featureTickets[tileNum] = frm.requestFeatures(dataSource,
                                                                 sourceID,
                                                                 listener);
                }
            }

            if (regKey != null) {
                featureTickets[tileNum].setFetchGroup(regKey);
            }
        }
    }

    void registerFeatureFetcher(Location loc, Object regKey) {
        Location[] tiles = getTiles();
        for (int t = 0; t < tiles.length; ++t) {
            if (LocationTools.overlaps(tiles[t], loc)) {
                registerFeatureFetcher(t, regKey);
            }
        }
    }

    void registerFeatureFetcher(Object regKey) {
        Location[] tiles = getTiles();
        for (int t = 0; t < tiles.length; ++t) {
            registerFeatureFetcher(t, regKey);
        }
    }

    protected FeatureHolder getFeatures() {
        if (allFeatures == null) {
            Location[] tiles = getTiles();
            if (tiles.length == 1) {
                allFeatures = new TileFeaturesWrapper(0);
            } else {
                try {
                    MergeFeatureHolder mfhAllFeatures = new MergeFeatureHolder();
                    for (int t = 0; t < tiles.length; ++t) {
                mfhAllFeatures.addFeatureHolder(new TileFeaturesWrapper(t));
                    }
                    mfhAllFeatures.addFeatureHolder(unrulyFeatures);
                    allFeatures = mfhAllFeatures;
                } catch (ChangeVetoException cve) {
                    throw new BioError(cve);
                }
            }
        }

        return allFeatures;
    }

    public Iterator features() {
        return getFeatures().features();
    }

    public boolean containsFeature(Feature f) {
      return getFeatures().containsFeature(f);
    }

    public FeatureHolder filter(FeatureFilter ff) {
        return filter(ff, !FilterUtils.areProperSubset(ff, FeatureFilter.top_level));
    }

    public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
        if (FilterUtils.areDisjoint(ff, getSchema())) {
            return FeatureHolder.EMPTY_FEATURE_HOLDER;
        } else {
            return getFeatures().filter(ff, recurse);
        }
    }

    public FeatureFilter getSchema() {
        FeatureFilter baseFilter = new FeatureFilter.And(
            new FeatureFilter.ByAnnotation(DASSequence.PROPERTY_ANNOTATIONSERVER,
                                                                           dataSource),
            getAllTypesFilter()
        );
        return new FeatureFilter.And(
                baseFilter,
                new FeatureFilter.Or(
                        new FeatureFilter.ByDescendant(baseFilter),
                        FeatureFilter.top_level
                )
        ) ;
    }

    public int countFeatures() {
        return getFeatures().countFeatures();
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
    // Listener for recieving the types document
    //

    private class DASTypesPopulator implements TypesListener {
        private Map types;

        public void startSegment() {
            types = new HashMap();
        }

        public void registerType(String type) {
            types.put(type, null);
        }

        public void registerType(String type, int count) {
            types.put(type, new Integer(count));
        }

        public void endSegment() {
            typesMap = types;
        }
    }

    //
    // Listener which is responsible for populating this FeatureSet
    //

    private class DASFeatureSetPopulator extends SeqIOAdapter {
        private SimpleFeatureHolder holder;
        private List featureStack = new ArrayList();
        private Feature stackTop = null;
        private int thisTile;
        private Location tileLocation = null;

        DASFeatureSetPopulator(int thisTile) {
            this.thisTile = thisTile;
            Location[] tiles = getTiles();
            if (tiles.length > 1) {
                this.tileLocation = tiles[thisTile];
            }
        }

        public void startSequence() {
            holder = new SimpleFeatureHolder();
        }

        public void endSequence() {
            tileFeatures[thisTile] = refSequence.getParentDB().getFeaturesCache().makeReference(holder);
            featureTickets[thisTile] = null;
        }

        public void startFeature(Feature.Template temp)
            throws ParseException
        {
            if (temp instanceof ComponentFeature.Template) {
                // I'm not convinced there's an easy, safe, way to say we don't
                // want these server side, so we'll elide them here instead.
                // We push a null onto the stack so that we don't get confused
                // over endFeature().

                featureStack.add(null);
            } else {
                try {
                    Feature f = null;
                    if (temp.annotation == Annotation.EMPTY_ANNOTATION) {
                        temp.annotation = new SmallAnnotation();
                    } else {
                        if (temp.annotation.containsProperty(XFFFeatureSetHandler.PROPERTY_XFF_ID)) {
                            temp.annotation.setProperty(DASSequence.PROPERTY_FEATUREID,
                                                        temp.annotation.getProperty(XFFFeatureSetHandler.PROPERTY_XFF_ID));
                        }
                    }
                    temp.annotation.setProperty(DASSequence.PROPERTY_ANNOTATIONSERVER, dataSource);

                    if (stackTop == null) {
                        f = ((RealizingFeatureHolder) refSequence).realizeFeature(refSequence, temp);

                        if (tileLocation == null || LocationTools.contains(tileLocation, f.getLocation())) {
                            holder.addFeature(f);
                        } else {
                            if (! unrulyFeatures.containsFeature(f)) {
                                unrulyFeatures.addFeature(f);
                            }
                        }
                    } else {
                        f = stackTop.createFeature(temp);
                    }

                    featureStack.add(f);
                    stackTop = f;
                } catch (Exception ex) {
                    ex.printStackTrace();
                    throw new ParseException(ex, "Couldn't realize feature in DAS");
                }
            }
        }

        public void addFeatureProperty(Object key, Object value)
            throws ParseException
        {
            if (stackTop == null) {
                // Feature we're skipping
                return;
            }

            try {
                if (key.equals(XFFFeatureSetHandler.PROPERTY_XFF_ID)) {
                    stackTop.getAnnotation().setProperty(DASSequence.PROPERTY_FEATUREID, value);
                } else {
                    Annotation ann = stackTop.getAnnotation();
                    if (ann.containsProperty(key)) {
                        Object o = ann.getProperty(key);
                        Collection col;
                        if (o instanceof Collection) {
                            col = (Collection) o;
                        } else {
                            col = new ArrayList();
                            col.add(o);
                            ann.setProperty(key, col);
                        }

                        col.add(value);
                    } else {
                        stackTop.getAnnotation().setProperty(key, value);
                    }
                }
            } catch (ChangeVetoException ex) {
                throw new ParseException(ex, "Couldn't set feature property");
            } catch (NullPointerException ex) {
                ex.printStackTrace();
            }
        }

        public void endFeature()
            throws ParseException
        {
            if (featureStack.size() < 1) {
                throw new BioError("Missmatched endFeature()");
            } else {
                featureStack.remove(featureStack.size() - 1);
                int pos = featureStack.size() - 1;
                stackTop = null;
                while (stackTop == null && pos >= 0) {
                    stackTop = (Feature) featureStack.get(pos--);
                }
            }
        }
    }

    private class TileFeaturesWrapper
      extends
        Unchangeable
      implements
        FeatureHolder
    {
        private int tileNum;

        TileFeaturesWrapper(int tileNum) {
            this.tileNum = tileNum;
        }

        protected FeatureHolder getFeatures() {
            if (tileFeatures[tileNum] != null) {
                FeatureHolder fh = (FeatureHolder) tileFeatures[tileNum].get();
                if (fh != null) {
                    return fh;
                }
            }

            registerFeatureFetcher(tileNum, null);
            try {
                featureTickets[tileNum].doFetch();
            } catch (Exception ex) {
                throw new BioRuntimeException(ex);
            }

            if (tileFeatures[tileNum] != null) {
                FeatureHolder fh = (FeatureHolder) tileFeatures[tileNum].get();
                if (fh != null) {
                    return fh;
                }
            }

            throw new BioRuntimeException("Feature fetch failed for now good reason...");
        }

        public int countFeatures() {
            return getFeatures().countFeatures();
        }

        public Iterator features() {
            return getFeatures().features();
        }

        public FeatureHolder filter(FeatureFilter ff, boolean recurse) {
            return getFeatures().filter(ff, recurse);
        }

    public FeatureHolder filter(FeatureFilter ff) {
        return getFeatures().filter(ff);
    }

        public Feature createFeature(Feature.Template templ)
            throws ChangeVetoException
        {
            throw new ChangeVetoException("NO");
        }

        public void removeFeature(Feature f)
            throws ChangeVetoException
        {
            throw new ChangeVetoException("NO");
        }

        public boolean containsFeature(Feature f) {
            return getFeatures().containsFeature(f);
        }

    public FeatureFilter getSchema() {
        FeatureFilter tileFilter = new FeatureFilter.ContainedByLocation(tiles[tileNum]);
        return new FeatureFilter.And(
                tileFilter,
                new FeatureFilter.OnlyDescendants(tileFilter)
        ) ;
    }
    }
}
