package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.io.PrintWriter;

import java.util.List;
import java.util.logging.Logger;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

public class CpgWalkerDomainFinderMethRange extends CpgWalkerDomainFinder {

	protected final int MIN_WEIGHTING_WIND_SIZE = 5000; // I admit, this is arbitrary and driven by CGI size
	
	protected double minMeth = 0.0;
	protected double maxMeth = 1.0;
//	protected boolean useWeighting = false;
	
	public CpgWalkerDomainFinderMethRange(CpgWalkerParams inWalkParams,
			String chr, PrintWriter pw, double inMinMeth, double inMaxMeth) {
		super(inWalkParams, chr, pw);
		this.minMeth = inMinMeth;
		this.maxMeth = inMaxMeth;
		
//		useWeighting = (inWalkParams.maxWindSize >= MIN_WEIGHTING_WIND_SIZE);
//		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe(
//				String.format("Use weighting if maxWindSize(%s) >= %d (useWeighting=%s)\n",inWalkParams.maxWindSize,MIN_WEIGHTING_WIND_SIZE,useWeighting));
	}

	@Override
	protected boolean windPasses(List<Cpg> inWindow) {

		boolean useWeighting = (CpgWalker.windLen(inWindow) >= this.MIN_WEIGHTING_WIND_SIZE);
		
		// We use the walker's built in meth summarizer
		double meth = this.methSummarizer.getValMean(useWeighting);
		
		if (this.walkParams.debug) System.err.printf("\tTesting window\t%s\n", CpgWalker.windStr(inWindow));

		boolean passes = ((meth>=this.minMeth) && (meth<=this.maxMeth));

		if (this.walkParams.debug && passes) System.err.printf("\t\tFound passing window\t%s\n", CpgWalker.windStr(inWindow));
//		if (this.walkParams.debug && passes) System.err.printf("\t\tFound passing window, meth =%.2f\tsize=%d (%d CpGs)\n",
//				meth, inWindow.get(inWindow.size()-1).chromPos-inWindow.get(0).chromPos, inWindow.size());
		

		return passes;
	}

	@Override
	protected double windScore(List<Cpg> inWindow) {
		boolean useWeighting = (CpgWalker.windLen(inWindow) >= this.MIN_WEIGHTING_WIND_SIZE);
		double meth = this.methSummarizer.getValMean(useWeighting);
		return meth;
	}

}
