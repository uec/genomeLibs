package edu.usc.epigenome.genomeLibs.GenomicRange;

import org.biojava.bio.seq.StrandedFeature;

public interface GenomicPositionScored {

	public StrandedFeature.Strand getStrand();
	public int getPos();
	public double getScore(StrandedFeature.Strand inStrand);


}
