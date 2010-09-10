package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.io.PrintWriter;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

public class CpgWalkerDomainFinderMethMeans extends CpgWalkerDomainFinder {

	protected final int MIN_WEIGHTING_WIND_SIZE = 5000; // I admit, this is arbitrary and driven by CGI size
	protected List<String> tables = null;
	protected boolean showCoords = false;
	
//	protected boolean useWeighting = false;
	
	public CpgWalkerDomainFinderMethMeans(CpgWalkerParams inWalkParams,
			String chr, PrintWriter pw, List<String> inTables) {
		super(inWalkParams, chr, pw, inTables.size());
		tables = inTables;
		myInit();

		
//		useWeighting = (inWalkParams.maxWindSize >= MIN_WEIGHTING_WIND_SIZE);
//		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe(
//				String.format("Use weighting if maxWindSize(%s) >= %d (useWeighting=%s)\n",inWalkParams.maxWindSize,MIN_WEIGHTING_WIND_SIZE,useWeighting));
	}

	
	
	protected void myInit() {

		String headers = null;
		try
		{
			headers = this.headerStr();
		}
		catch (Exception e)
		{
			System.err.printf("CpgWalker::init() exception: %s\n",e);
			e.printStackTrace();
		}
		
		if (headers!= null)
		{
			pw.println(headers);
		}
	}



	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#headerStr()
	 */
	@Override
	public String headerStr() throws Exception {
		String out = null;
		if (this.tables != null)
		{
			List<String> list = new ArrayList<String>(this.tables);
			list.add(0, "nCpgs");
			if (this.showCoords)
			{
				list.add(0, "end");
				list.add(0, "start");
				list.add(0, "chr");
			}
			ListUtils.setDelim(",");
			out= ListUtils.excelLine(list);
		}
		return out;
	}



	@Override
	protected boolean windPasses(List<Cpg[]> inWindow) {

		boolean useWeighting = (CpgWalker.windLen(inWindow, true) >= this.MIN_WEIGHTING_WIND_SIZE);
		
		// We use the walker's built in meth summarizer
		double meths[] = new double[this.numTables()];
		int start = CpgWalker.windStart(inWindow, false);
		int end = CpgWalker.windEnd(inWindow, false);
		int nCpgs = inWindow.size();
		String chrom = this.getCurChr();
		if (this.showCoords) this.pw.printf("%s,%d,%d,", chrom,start,end);
		this.pw.printf("%d,",nCpgs);
		
		for (int t = 0; t < this.numTables(); t++)
		{
			//System.err.printf("Adding to ALL printwriter, %s\n",CpgWalker.windStrMulti(inWindow));
			meths[t] = this.methSummarizer.get(t).getValMean(useWeighting);
			this.pw.printf("%s%.4f", ((t>0)?",":""), meths[t]);
		}
		this.pw.println();

		return false;
	}

	
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerDomainFinder#outputCurWind(edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange)
	 * 
	 * Don't output
	 */
	@Override
	protected void outputCurWind(GenomicRange gr) {
	}

	@Override
	protected double windScore(List<Cpg[]> inWindow) {
		boolean useWeighting = (CpgWalker.windLen(inWindow, true) >= this.MIN_WEIGHTING_WIND_SIZE);
		double meth = this.methSummarizer.get(0).getValMean(useWeighting);
		return meth;
	}

}
