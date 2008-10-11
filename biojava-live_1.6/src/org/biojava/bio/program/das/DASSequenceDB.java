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

import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

import javax.xml.parsers.DocumentBuilder;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.MergeFeatureHolder;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.db.SequenceDBLite;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.Unchangeable;
import org.biojava.utils.cache.Cache;
import org.biojava.utils.cache.FixedSizeCache;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

/**
 * <p>
 * Collection of sequences retrieved from the DAS network.
 * </p>
 *
 * <p>
 * The DAS-specific parts of this API are still subject
 * to change.
 * </p>
 *
 * @author Thomas Down
 * @since 1.1
 */

public class DASSequenceDB
  extends
    Unchangeable
  implements
    SequenceDB
{
    private URL dataSourceURL;
    private Map sequences;
    private Cache symbolsCache;
    private FixedSizeCache featuresCache;
    private Set rootIDs;
    private FeatureRequestManager frm;

    {
        sequences = new HashMap();
        symbolsCache = new FixedSizeCache(20);
        featuresCache = new FixedSizeCache(50);
    }

    Cache getSymbolsCache() {
        return symbolsCache;
    }

    /**
     * @throws BioException if the capacity can't be reached.
     */
    void ensureFeaturesCacheCapacity(int min) throws BioException {
      //  if(min > MAX_CAPACITY) {
//          throw new BioException( "Capacity of (" + MAX_CAPACITY +
//                                  " exceeded by " + min);
//        }
        if (featuresCache.getLimit() < min) {
            // System.err.println("Setting cache limit up to " + min);
            featuresCache.setLimit(min);
        }
    }

    Cache getFeaturesCache() {
        return featuresCache;
    }

    FeatureRequestManager getFeatureRequestManager() {
        if (frm == null) {
            frm = new FeatureRequestManager(this);
        }

        return frm;
    }

    DASSequenceDB() {
        // Constructor for dummy objects.  Ugh.
    }


  public FeatureHolder filter(FeatureFilter ff) {
      MergeFeatureHolder results = new MergeFeatureHolder();
      try {
          for (SequenceIterator si = sequenceIterator(); si.hasNext(); ) {
              Sequence seq = si.nextSequence();
              FeatureHolder fh = seq.filter(ff);
              if (fh != FeatureHolder.EMPTY_FEATURE_HOLDER) {
                  results.addFeatureHolder(fh);
              }
          }
      } catch (BioException ex) {
          throw new BioRuntimeException(ex);
      } catch (ChangeVetoException cve) {
          throw new BioError("Assertion failed: couldn't modify newly created MergeFeatureHolder", cve);
      }
      return results;
  }

    public DASSequenceDB(URL dataSourceURL)
        throws BioException
    {
        String s = dataSourceURL.toString();
        if (! (s.endsWith("/"))) {
            try {
                dataSourceURL = new URL(s + "/");
            } catch (MalformedURLException ex) {
                throw new BioException("Assertion failure: trivial URI manipulation failed",ex);
            }
        }
        this.dataSourceURL = dataSourceURL;
    }

    DASSequence _getSequence(String id)
        throws BioException, IllegalIDException
    {
        return _getSequence(id, Collections.singleton(dataSourceURL));
    }

    DASSequence _getSequence(String id, Set annoURLs)
        throws BioException, IllegalIDException
    {
        DASSequence seq = (DASSequence) sequences.get(id);
        if (seq == null) {
            seq = new DASSequence(this, dataSourceURL, id, annoURLs);
            sequences.put(id, seq);
        }
        return seq;
    }

    /**
     * Return a SequenceDB exposing /all/ the entry points
     * in this DAS datasource.
     *
     */

    private SequenceDBLite allEntryPoints;

    public SequenceDBLite allEntryPointsDB() {
        if (allEntryPoints == null) {
            allEntryPoints = new AllEntryPoints();
        }

        return allEntryPoints;
    }

    private class AllEntryPoints
      extends
        Unchangeable
      implements
        SequenceDBLite
    {
        public Sequence getSequence(String id)
            throws BioException, IllegalIDException
        {
            return _getSequence(id);
        }

        public void addSequence(Sequence seq)
            throws ChangeVetoException
        {
            throw new ChangeVetoException("No way we're adding sequences to DAS");
        }

        public void removeSequence(String id)
            throws ChangeVetoException
        {
            throw new ChangeVetoException("No way we're removing sequences from DAS");
        }

        public String getName() {
            return "All sequences in " + dataSourceURL.toString();
        }
  }


    /**
     * Return the URL of the reference server for this database.
     */

    public URL getURL() {
        return dataSourceURL;
    }

    public String getName() {
        return dataSourceURL.toString();
    }

    public Sequence getSequence(String id)
        throws BioException, IllegalIDException
    {
        if (! (ids().contains(id))) {
            throw new IllegalIDException("Database does not contain " + id + " as a top-level sequence");
        }
        return _getSequence(id);
    }

    public Set ids() {
        if (rootIDs == null) {
            try {
                DAS.startedActivity(this);

                Set ids = new HashSet();

                URL epURL = new URL(dataSourceURL, "entry_points");
                HttpURLConnection huc = (HttpURLConnection) epURL.openConnection();
                try {
                    huc.connect();
                } catch (Exception e) {
                    throw new BioException("Can't connect to " + epURL, e);
                }
                int status = DASSequenceDB.tolerantIntHeader(huc, "X-DAS-Status");
                if (status == 0)
                    throw new BioException("Not a DAS server: " + dataSourceURL + " Query: " + epURL);
                else if (status != 200)
                    throw new BioException("DAS error (status code = " + status +
                                           ") connecting to " + dataSourceURL + " with query " + epURL);


                InputSource is = new InputSource(huc.getInputStream());
                is.setSystemId(epURL.toString());
                DocumentBuilder parser = DASSequence.nonvalidatingParser();
                Element el = parser.parse(is).getDocumentElement();

                NodeList segl = el.getElementsByTagName("SEGMENT");
                for (int i = 0; i < segl.getLength(); ++i) {
                    el = (Element) segl.item(i);
                    String id = el.getAttribute("id");
                    ids.add(id);
                }

                rootIDs = Collections.unmodifiableSet(ids);
            } catch (SAXException ex) {
                throw new BioRuntimeException("Exception parsing DAS XML",ex);
            } catch (IOException ex) {
                throw new BioRuntimeException("Error connecting to DAS server",ex);
            } catch (NumberFormatException ex) {
                throw new BioRuntimeException("Error parsing number",ex);
            } catch (BioException ex) {
                throw new BioRuntimeException(ex);
            } finally {
                DAS.completedActivity(this);
            }
        }

        return rootIDs;
    }

    public void addSequence(Sequence seq)
        throws ChangeVetoException
    {
        throw new ChangeVetoException("No way we're adding sequences to DAS");
    }

    public void removeSequence(String id)
        throws ChangeVetoException
    {
        throw new ChangeVetoException("No way we're removing sequences from DAS");
    }

    public SequenceIterator sequenceIterator()
    {
        return new SequenceIterator() {
            private Iterator i = ids().iterator();

            public boolean hasNext() {
                return i.hasNext();
            }

            public Sequence nextSequence()
                throws BioException
            {
                return getSequence((String) i.next());
            }
        } ;
    }

    static int tolerantIntHeader(HttpURLConnection huc, String name)
    {
        try {
            String header = huc.getHeaderField(name);
            if (header == null) {
                return 0;
            }

            String firstToken = new StringTokenizer(header).nextToken();
            return Integer.parseInt(firstToken);
        } catch (NumberFormatException ex) {
            return 0;
        }
    }
}
