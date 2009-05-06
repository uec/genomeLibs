/**
 * 
 */
package edu.usc.epigenome.genomeLibs.GenomicRange;

import java.util.Comparator;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.BiojavaUtils;
import edu.usc.epigenome.genomeLibs.ChromStringComparator;

/**
 * @author benb
 * 
 * Correctly sorts genomic ranges where chromosomes have the form
 * c(hr)xx where "xx" is a number.  So order will be chr1, chr2,
 * chr3 (or c1, c2, c3), rather than the native (and much faster 
 * 
 */

public class GenomicRangeComparator implements Comparator<GenomicRange> {

	protected boolean fastSort = true;

	/**
	 * @param inFastSort If true, the sort is fast but chromsomes will be
	 * sorted chr1, chr10, chr11, ..., chr2, chr20, ...
	 */
	public GenomicRangeComparator(boolean inFastSort) {
		super();
		this.fastSort = inFastSort;
	}

		
	public int compare(GenomicRange t, GenomicRange o) {

		return compareHelper(t, o, this.fastSort);
	}

	public static int compareHelper(GenomicRange t, GenomicRange o, boolean fast) {
		
		int strCmp = (fast) ? t.getChrom().compareTo(o.getChrom()) :  
			(new ChromStringComparator()).compare(t.getChrom(), o.getChrom());

		
		if (strCmp != 0) return strCmp;

		// Do the end first (this only matters if one completely 
		// subsumes the other, in which case the internal one will
		// be returned first.
		if (t.getEnd() > o.getEnd()) return 1; 
		if (t.getEnd() < o.getEnd()) return -1;
		
		if (t.getStart() > o.getStart()) return 1; 
		if (t.getStart() < o.getStart()) return -1;
		
		if (BiojavaUtils.strandToInt(t.getStrand()) > BiojavaUtils.strandToInt(o.getStrand())) return 1; 
		if (BiojavaUtils.strandToInt(t.getStrand()) < BiojavaUtils.strandToInt(o.getStrand())) return -1; 

		return 0;
	}	
		
}
