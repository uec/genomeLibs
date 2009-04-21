package edu.usc.epigenome.genomeLibs.GenomicRange;

import org.biojava.bio.seq.StrandedFeature;

public interface GenomicPositionScored {

	public StrandedFeature.Strand getStrand();
	public String getChr();
	public int getPos();
	public double getSummaryScore();
	public double getStrandedScore(StrandedFeature.Strand inStrand);

}
