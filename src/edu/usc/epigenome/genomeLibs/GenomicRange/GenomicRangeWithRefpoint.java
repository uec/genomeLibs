package edu.usc.epigenome.genomeLibs.GenomicRange;

import org.biojava.bio.seq.StrandedFeature.Strand;

public class GenomicRangeWithRefpoint extends GenomicRange {
	
	int refPoint = -1;

	public GenomicRangeWithRefpoint(String inChrom, int inStart, int inEnd, int inRefPoint) {
		super(inChrom, inStart, inEnd);
		this.setRefPoint(inRefPoint);
	}

	public GenomicRangeWithRefpoint(String inChrom, int inStart, int inEnd) {
		super(inChrom, inStart, inEnd);
	}

	public GenomicRangeWithRefpoint(String inChrom, int inStart, int inEnd,
			Strand inStrand) {
		super(inChrom, inStart, inEnd, inStrand);
	}

	public int getRefPoint() {
		return (refPoint == -1) ? this.start : refPoint;
	}

	public void setRefPoint(int refPoint) {
		this.refPoint = refPoint;
	}
	
	
	

}
