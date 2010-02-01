package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.io.PrintWriter;

public class CpgWalkerDomainFinderMethRange extends CpgWalkerDomainFinder {

	protected double minMeth = 0.0;
	protected double maxMeth = 1.0;
	
	public CpgWalkerDomainFinderMethRange(CpgWalkerParams inWalkParams,
			String chr, PrintWriter pw, double inMinMeth, double inMaxMeth) {
		super(inWalkParams, chr, pw);
		this.minMeth = inMinMeth;
		this.maxMeth = inMaxMeth;
		
	}

	@Override
	protected boolean windPasses() {
		
		// We use the walker's built in meth summarizer
		double meth = this.methSummarizer.getValMean(true);
		return ((meth>=this.minMeth) && (meth<=this.maxMeth));
	}

	@Override
	protected double windScore() {
		double meth = this.methSummarizer.getValMean();
		return meth;
	}

}
