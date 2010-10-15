package edu.usc.epigenome.scripts;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import com.googlecode.charts4j.AxisLabels;
import com.googlecode.charts4j.AxisLabelsFactory;
import com.googlecode.charts4j.BarChart;
import com.googlecode.charts4j.BarChartPlot;
import com.googlecode.charts4j.Color;
import com.googlecode.charts4j.Data;
import com.googlecode.charts4j.DataEncoding;
import com.googlecode.charts4j.DataUtil;
import com.googlecode.charts4j.Fills;
import com.googlecode.charts4j.GCharts;
import com.googlecode.charts4j.Plot;
import com.googlecode.charts4j.Plots;
import com.googlecode.charts4j.VennDiagram;
import com.sun.org.apache.xalan.internal.xsltc.compiler.Pattern;
import com.sun.tools.javac.code.Attribute.Array;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.FeatDb.FeatIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgCoverageSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgDeaminationSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgDensitySummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgNonconversionSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgSummarizer;



public class MethylDbToFeatPairCounts {

	private static final String C_USAGE = "Use: MethylDbToFeatPairCounts featType1 featType2 .. ";
	
	public int STEP = (int)1E6;
	public int MAXCOORD = (int)2.8E8;
	
	// receives other command line parameters than options
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-threeWay",usage="Do feature triplets")
    protected boolean threeWay = false;
    @Option(name="-fourWay",usage="Do feature quadruplets")
    protected boolean fourWay = false;
    @Option(name="-cpgCounts",usage="Do counts based on raw number of Cpgs (rather than base pairs)")
    protected boolean cpgCounts = false;
    @Option(name="-outPrefix",usage="outputFilename")
    protected String outPrefix = "featPairs";
    @Option(name="-masterFeat",usage="If a master feat is specified, we only do pairwise against that")
    protected String masterFeat = null;
    @Option(name="-methylDbPrefix",usage="use this table to get CpGs")
    protected String methylDbPrefix = "methylCGsRich_normalM030510_";
    @Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
    protected int minCTreads = 0;  // For this, we just care about whole genome stats
    @Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
    		"This sets a maximum number of observed As on the opposite strand (default 0.1)")
    protected double maxOppStrandAfrac = 0.1;
    
	@Argument
	private List<String> featTypes = new ArrayList<String>();
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToFeatPairCounts().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception
	{
		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);

			if(featTypes.size() < 1 ) {
				System.err.println(C_USAGE);
				parser.printUsage(System.err);
				System.exit(1);
			}

			if ((masterFeat != null) && (this.threeWay || this.fourWay))
			{
				System.err.println("Can not do three or four way if masterFeat is specified");
				System.err.println(C_USAGE);
				parser.printUsage(System.err);
				System.exit(1);
			}

		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			System.err.println(C_USAGE);
			// print the list of available options
			parser.printUsage(System.err);
			System.err.println();
			return;
		}
		
		int nFeats = featTypes.size();
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.SEVERE);
		String countsSec = (this.cpgCounts) ? "cpgCounts" : "bpCounts"; 
		String vennFn = String.format("Venns%s.%s.htm", this.outPrefix, countsSec);
		PrintWriter vennPw = new PrintWriter(new FileOutputStream(vennFn));
		
		String csvFn = String.format("Venns%s.%s.csv", this.outPrefix, countsSec);
		PrintWriter csvPw = new PrintWriter(new FileOutputStream(csvFn));
		
		
		// Go through individual feats.
		// counts can be larger than 32-bit Integer limit
		Map<String,Long> savedCounts = new HashMap<String,Long>(nFeats);
		List<String> singles = new ArrayList<String>(nFeats+1);
		if (this.masterFeat!=null) singles.add(this.masterFeat);
		singles.addAll(featTypes);
		for (String featType : singles)
		{
			List<String> feats = new ArrayList<String>(1);
			feats.add(featType);
			long count = calculateAndPrintCounts(feats, savedCounts); // , vennPw, csvPw);
		}

		// Then feat pairs
		if (this.masterFeat != null)
		{
			for (int i = 0; i < nFeats; i++)
			{
				List<String> feats = Arrays.asList(this.masterFeat, featTypes.get(i));
				long count = calculateAndPrintCounts(feats, savedCounts, vennPw, csvPw);
			}
		}
		else
		{
			if (nFeats>1)
			{
				for (int i = 0; i < nFeats; i++)
				{
					for (int j = i+1; j < nFeats; j++)
					{
						List<String> feats = new ArrayList<String>(2);
						feats.add(featTypes.get(i));
						feats.add(featTypes.get(j));
						long count = calculateAndPrintCounts(feats, savedCounts, vennPw, csvPw);
					}
				}
			}

			// Triplets. I know there's a fancy recursive way to do this, but 
			// this produces a nice ordering for the time being.
			if (threeWay && (nFeats>2))
			{
				for (int i = 0; i < nFeats; i++)
				{
					for (int j = i+1; j < nFeats; j++)
					{
						for (int k = j+1; k < nFeats; k++)
						{
							List<String> feats = new ArrayList<String>(3);
							feats.add(featTypes.get(i));
							feats.add(featTypes.get(j));
							feats.add(featTypes.get(k));
							calculateAndPrintCounts(feats, savedCounts, vennPw, csvPw);
						}
					}
				}
			}

			// Don't have a way to draw quadruplets, but we can output them
			if (fourWay && (nFeats>3))
			{
				for (int i = 0; i < nFeats; i++)
				{
					for (int j = i+1; j < nFeats; j++)
					{
						for (int k = j+1; k < nFeats; k++)
						{
							for (int l = k+1; l < nFeats; l++)
							{
								List<String> feats = new ArrayList<String>(3);
								feats.add(featTypes.get(i));
								feats.add(featTypes.get(j));
								feats.add(featTypes.get(k));
								feats.add(featTypes.get(l));
								// Put null since we can't do the venns
								calculateAndPrintCounts(feats, savedCounts, null,csvPw);
							}
						}
					}
				}
			}
		}
		
		vennPw.close();
		csvPw.close();
	}

	protected long calculateAndPrintCounts(List<String> featTypes, Map<String,Long> savedCounts)
	throws Exception
	{
		return calculateAndPrintCounts(featTypes, savedCounts, null, null);
	}
	
	protected long calculateAndPrintCounts(List<String> featTypes, Map<String,Long> savedCounts, 
			PrintWriter vennPw, PrintWriter csvPw) 
	throws Exception
	{
		StringBuffer sb = new StringBuffer(10000);
		int nFeats = featTypes.size();
		
		System.err.printf("About to calculateAndPrintCounts for: %s\n", ListUtils.excelLine(featTypes));
		
		// And the querier params
		MethylDbQuerier params = new MethylDbQuerier();
		params.methylTablePrefix = this.methylDbPrefix;
		params.setMinCTreads(this.minCTreads);
		params.setUseNonconversionFilter(!this.noNonconvFilter);
		params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);

		// Feat filts
		params.addFeatFilters(featTypes);
		
		ListUtils.setDelim(",");
		String featStr = ListUtils.excelLine(featTypes);
		
		// Count can be over the INT 32-bit limit
		long count = 0;
		
		for (String chr :  Arrays.asList("chrX")) // MethylDbUtils.CHROMS) //  Arrays.asList("chr21")) //
		{		
			
			System.err.printf("Feats=%s\tchr=%s\n",featStr, chr);
			// Iterator uses DB connection and can use a ton of memory because
			// it loads all rows at once.  This stuff should really be added to iterator
			// class, but until it is , just iterate here over the chromosome
			int MINCOORD = 0;
			MAXCOORD = GoldAssembly.chromLengthStatic(chr, "hg18");
			
//			MINCOORD = 8000000;
//			MAXCOORD = 10000000;

			for (int c = MINCOORD; c < MAXCOORD; c += STEP)
			{
				int last = Math.min(c+STEP-1, MAXCOORD);
				params.clearRangeFilters();
				params.addRangeFilter(chr, c, last);
//				params.addRangeFilter(chr);
			
				CpgIterator cpgit = new CpgIterator(params);
				
				// Weighted or unweighted
				long counti;
				if (this.cpgCounts)
				{
					int numRows = cpgit.getCurNumRows();
					counti = numRows;
				}
				else
				{
					counti = 0;
					while (cpgit.hasNext())
					{
						Cpg cpg = cpgit.next();
						counti += (int)cpg.getCpgWeight();
						//System.err.printf("Weight=%d\n", (int)cpg.getCpgWeight());
					}
				}
				
				count += counti;
			}
		}
		

		if ((nFeats > 1) && (vennPw != null))
		{
			vennPw.printf("<H3>%s</H3>\n", featStr);
			if (count==0)
			{
				vennPw.printf("<H4>0 base pairs in common</H4>\n");
			}
			else
			{

			
			// If you make them too long
			String a = featTypes.get(0) + "(" + savedCounts.get(featTypes.get(0)) + ")"; 
			String b = featTypes.get(1) + "(" + savedCounts.get(featTypes.get(1)) + ")";
			String c = (nFeats>2) ? (featTypes.get(2) + "(" + savedCounts.get(featTypes.get(2)) + ")") : "";

			double an = savedCounts.get(featTypes.get(0));
			double bn = savedCounts.get(featTypes.get(1));
			double cn = (nFeats>2) ? savedCounts.get(featTypes.get(2)) : 0;
			double abn = (nFeats==2) ? count : savedCounts.get(featTypes.get(0) + "+" + featTypes.get(1));
			double acn = (nFeats==2) ? 0 : savedCounts.get(featTypes.get(0) + "+" + featTypes.get(2));
			double bcn = (nFeats==2) ? 0 : savedCounts.get(featTypes.get(1) + "+" + featTypes.get(2));
			double abcn = (nFeats==2) ? 0 : count;
			double[] allVals = {an, bn, cn, abn, acn, bcn, abcn};
			ListUtils.setDelim(", ");
			System.err.println("vals: " + ListUtils.excelLine(allVals));
			double divBy = MatUtils.nanMax(allVals) / 99.99; // Sometimes there's an overflow to 100.00001 when i set this to 100.0, and charts4j crashes

			
//			VennDiagram chart = GCharts.newVennDiagram(100, 80, 60, 30, 30, 30, 10);
			VennDiagram chart = GCharts.newVennDiagram(an/divBy, bn/divBy, cn/divBy, 
					abn/divBy, acn/divBy, bcn/divBy, abcn/divBy);
			chart.setDataEncoding(DataEncoding.TEXT);
			chart.setTitle(featStr, Color.WHITE, 16);
			chart.setSize(600, 200);
			chart.setCircleLegends(a,b,c);
			chart.setCircleColors(Color.GREEN, Color.RED, Color.BLUE);
			//chart.setBackgroundFill(Fills.newSolidFill(Color.BLACK));
			String url = chart.toURLString();
			// EXAMPLE CODE END. Use this url string in your web or
			// Internet application.

			vennPw.printf("<IMG SRC=\"%s\">\n", url);
		}
		}


		// CSV Output
		if (csvPw != null)
		{
			sb.append(String.format("%s",featStr));
			for (int i = 0; i < nFeats; i++)
			{
				Long singleCount = savedCounts.get(featTypes.get(i));
				if (singleCount != null)
				{
					Long remainder = singleCount-count;
					double sens = (double)count/(double)singleCount;
					sb.append(String.format(",%d,%.4f", count, sens)); 
					sb.append(String.format(",%d,%.4f", remainder, 1-sens)); 
				}
			}

			sb.append("\n");
			System.err.println("Setting " + featStr + ": " + count);
			csvPw.append(sb);
		}

		savedCounts.put(featStr, count);

		return count;
	}
	


	
}

