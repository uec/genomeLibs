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
import java.util.Collections;
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.distributed.DistDataSource;

/**
 * <p>
 * View of DAS data suitable for integration via the
 * meta-DAS system.
 * </p>
 *
 * <p>
 * This class represents an alternative view on DAS data, designed
 * to be used with the new MetaDAS data integration framework, rather
 * than including built-in integration code.  It shares all the
 * network code, and much other code, with the older
 * <code>DASSequence</code> API.  This new API is currently quite
 * experimental, so stick to <code>DASSequence</code> for now...
 * </p>
 *
 * @author Thomas Down
 * @since 1.2 [MetaDAS] 
 */

class DASAnnotationDistDataSource implements DistDataSource {
    private URL datasource;
    private DASSequenceDB dummyDB = new DASSequenceDB();

    public URL getURL() {
	return datasource;
    }

    public DASAnnotationDistDataSource(URL url) 
        throws BioException
    {
	this.datasource = url;
    }

    public boolean hasSequence(String id) throws BioException {
	return false;
    }

    public boolean hasFeatures(String id) throws BioException {
	try {
	    getFeatures(id, FeatureFilter.all, false);
	    return true;
	} catch (Exception ex) {
	    return false;
	}
    }

    public FeatureHolder getFeatures(FeatureFilter ff) throws BioException {
	throw new BioException();
    }

    public FeatureHolder getFeatures(String id, FeatureFilter ff, boolean recurse) throws BioException {
	FeatureHolder fh;
	try {
	    fh = new RawAnnotationSequence(dummyDB, datasource, id);
	} catch (IllegalIDException ex) {
	    return FeatureHolder.EMPTY_FEATURE_HOLDER;
	}
	
	if (recurse == false && FilterUtils.areProperSubset(FeatureFilter.all, ff)) {
	    return fh;
	} else {
	    return fh.filter(ff, recurse);
	}
    }

    public Sequence getSequence(String id) throws BioException {
	throw new BioException("Erk");
    }

    public Set ids(boolean topLevel) throws BioException {
	return Collections.EMPTY_SET;
    }

    public boolean equals(Object o) {
	if (o instanceof DASAnnotationDistDataSource) {
	    return ((DASAnnotationDistDataSource) o).getURL().equals(getURL());
	}
	return false;
    }

    public int hashCode() {
	return getURL().hashCode() + 5;
    }
}
