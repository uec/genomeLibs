package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.io.PrintWriter;

import java.util.List;
import java.util.logging.Logger;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

public class CpgWalkerDomainFinderMethDiffs extends CpgWalkerDomainFinder {

	protected final int MIN_WEIGHTING_WIND_SIZE = 5000; // I admit, this is arbitrary and driven by CGI size
	
	protected double tableLowMethMaxMeth = 1.01;
	protected double tableHighMethMinMeth = -0.01;
	protected double minMethDiff = 0.0;
	protected int lowMethIndex = 0;
	protected int highMethIndex = 1;
//	protected boolean useWeighting = false;
	
	public CpgWalkerDomainFinderMethDiffs(CpgWalkerParams inWalkParams,
			String chr, PrintWriter pw, double inTableLowMethMaxMeth, double inTableHighMethMinMeth,
			double inMinMethDiff, int inLowMethIndex, int inHighMethIndex) {
		super(inWalkParams, chr, pw,2);
		
		if (this.numTables()!=2)
		{
			System.err.printf("CpgWalkerDomainFinderMethDiffs instantiated with %d tables, when it should be exactly two. Exiting\n", this.numTables());
			System.exit(1);
		}
		
		this.lowMethIndex = inLowMethIndex;
		this.highMethIndex = inHighMethIndex;
		this.tableLowMethMaxMeth = inTableLowMethMaxMeth;
		this.tableHighMethMinMeth = inTableHighMethMinMeth;
		this.minMethDiff = inMinMethDiff;
		
//		useWeighting = (inWalkParams.maxWindSize >= MIN_WEIGHTING_WIND_SIZE);
//		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe(
//				String.format("Use weighting if maxWindSize(%s) >= %d (useWeighting=%s)\n",inWalkParams.maxWindSize,MIN_WEIGHTING_WIND_SIZE,useWeighting));
	}

	@Override
	protected boolean windPasses(List<Cpg[]> inWindow) {

		boolean useWeighting = (CpgWalker.windLen(inWindow, true) >= this.MIN_WEIGHTING_WIND_SIZE);
		
		// We use the walker's built in meth summarizer
		double methlow = this.methSummarizer.get(this.lowMethIndex).getValMean(useWeighting);
		double methhigh = this.methSummarizer.get(this.highMethIndex).getValMean(useWeighting);
		double diff = methhigh - methlow;
		
//		if (this.walkParams.debug) System.err.printf("\tTesting window\t%s\n\tmethlow=%.2f\tmethhigh=%.2f\n", 
//				CpgWalker.windStr(inWindow,false),methlow,methhigh);

		boolean passes = 
			(methlow <= this.tableLowMethMaxMeth) && 
			(methhigh >= this.tableHighMethMinMeth) &&
			(diff >= this.minMethDiff);

		if (this.walkParams.debug && passes) System.err.printf("\t\tFound passing window\t%s\n", CpgWalker.windStr(inWindow,false));
//		if (this.walkParams.debug && passes) System.err.printf("\t\tFound passing window, meth =%.2f\tsize=%d (%d CpGs)\n",
//				meth, inWindow.get(inWindow.size()-1).chromPos-inWindow.get(0).chromPos, inWindow.size());
		

		return passes;
	}

	@Override
	protected double windScore(List<Cpg[]> inWindow) {
		boolean useWeighting = (CpgWalker.windLen(inWindow, true) >= this.MIN_WEIGHTING_WIND_SIZE);
		double meth = this.methSummarizer.get(0).getValMean(useWeighting);
		return meth;
	}

}
