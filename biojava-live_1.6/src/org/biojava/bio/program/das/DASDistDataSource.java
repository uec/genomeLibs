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
import java.util.Set;

import org.biojava.bio.BioException;
import org.biojava.bio.seq.FeatureFilter;
import org.biojava.bio.seq.FeatureHolder;
import org.biojava.bio.seq.FilterUtils;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.db.IllegalIDException;
import org.biojava.bio.seq.distributed.DistDataSource;

/**
 * View of DAS data suitable for integration via the
 * meta-DAS system.
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

class DASDistDataSource implements DistDataSource {
    private DASSequenceDB db;

    public URL getURL() {
	return db.getURL();
    }

    public DASDistDataSource(URL url) 
        throws BioException
    {
	this.db = new DASSequenceDB(url);
    }

    public boolean hasSequence(String id) throws BioException {
	try {
	    Sequence seq = db.allEntryPointsDB().getSequence(id);
	    int size = seq.length();
	    size = size+1-1; // trick to prevent unused field warning
	} catch (IllegalIDException ex) {
	    return false;
	}
	return true;
    }

    public boolean hasFeatures(String id) throws BioException {
	return hasSequence(id);
    }

    public FeatureHolder getFeatures(FeatureFilter ff) throws BioException {
	throw new BioException();
    }

    public FeatureHolder getFeatures(String id, FeatureFilter ff, boolean recurse) throws BioException {
	FeatureHolder fh;
	try {
	    fh = db.allEntryPointsDB().getSequence(id);
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
	return db.allEntryPointsDB().getSequence(id);
    }

    public Set ids(boolean topLevel) throws BioException {
	return db.ids();
    }
}
