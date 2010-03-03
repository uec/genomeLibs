package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.io.PrintWriter;

import java.util.List;
import java.util.logging.Logger;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

public class CpgWalkerDomainFinderMethRange extends CpgWalkerDomainFinder {

	protected final int MIN_WEIGHTING_WIND_SIZE = 20000;
	
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
	protected boolean windPasses(List<Cpg> inWindow) {
		
		// We use the walker's built in meth summarizer
		double meth = this.methSummarizer.getValMean(this.useWeighting);
		
		//System.err.printf("\tChecking window, meth =%.2f\n",meth);
		
		return ((meth>=this.minMeth) && (meth<=this.maxMeth));
	}

	@Override
	protected double windScore(List<Cpg> inWindow) {
		double meth = this.methSummarizer.getValMean(this.useWeighting);
		return meth;
	}

}
