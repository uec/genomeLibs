package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

public class CpgWalkerDomainFinderMethRange extends CpgWalkerDomainFinder {

	protected double minMeth = 0.0;
	protected double maxMeth = 1.0;
	
	public CpgWalkerDomainFinderMethRange(CpgWalkerParams inWalkParams,
			String chr, double inMinMeth, double inMaxMeth) {
		super(inWalkParams, chr);
		this.minMeth = inMinMeth;
		this.maxMeth = inMaxMeth;
	}

	@Override
	protected boolean windPasses() {
		
		// We use the walker's built in meth summarizer
		double meth = this.methSummarizer.getValMean();
		return ((meth>=this.minMeth) && (meth<=this.maxMeth));
	}

	@Override
	protected double windScore() {
		double meth = this.methSummarizer.getValMean();
		return meth;
	}

}
