package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.io.PrintWriter;
import java.util.logging.Logger;

public class CpgWalkerDomainFinderMethRange extends CpgWalkerDomainFinder {

	protected final int MIN_WEIGHTING_WIND_SIZE = 10000;
	
	protected double minMeth = 0.0;
	protected double maxMeth = 1.0;
	protected boolean useWeighting = false;
	
	public CpgWalkerDomainFinderMethRange(CpgWalkerParams inWalkParams,
			String chr, PrintWriter pw, double inMinMeth, double inMaxMeth) {
		super(inWalkParams, chr, pw);
		this.minMeth = inMinMeth;
		this.maxMeth = inMaxMeth;
		
		useWeighting = (inWalkParams.maxWindSize >= MIN_WEIGHTING_WIND_SIZE);
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe(
				String.format("Use weighting if maxWindSize(%s) >= %d (useWeighting=%s)\n",inWalkParams.maxWindSize,MIN_WEIGHTING_WIND_SIZE,useWeighting));
	}

	@Override
	protected boolean windPasses() {
		
		// We use the walker's built in meth summarizer
		double meth = this.methSummarizer.getValMean(this.useWeighting);
		
	//	System.err.printf("Checking window, meth =%.2f\n",meth);
		

		return ((meth>=this.minMeth) && (meth<=this.maxMeth));
	}

	@Override
	protected double windScore() {
		double meth = this.methSummarizer.getValMean(this.useWeighting);
		return meth;
	}

}
