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
    @Option(name="-methylDbPrefix",usage="use this table to get CpGs")
    protected String methylDbPrefix = "methylCGsRich_normal010310_";
    @Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
    protected int minCTreads = 2;
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
		String vennFn = String.format("Venns%s.htm", featTypes.get(0));
		PrintWriter vennPw = new PrintWriter(new FileOutputStream(vennFn));
		
		// Go through individual feats.
		Map<String,Integer> savedCounts = new HashMap<String,Integer>(nFeats);
		for (String featType : featTypes)
		{
			List<String> feats = new ArrayList<String>(1);
			feats.add(featType);
			int count = printCounts(feats, savedCounts, vennPw);
		}

		// Then feat pairs
		if (nFeats>1)
		{
			for (int i = 0; i < nFeats; i++)
			{
				for (int j = i+1; j < nFeats; j++)
				{
					List<String> feats = new ArrayList<String>(2);
					feats.add(featTypes.get(i));
					feats.add(featTypes.get(j));
					int count = printCounts(feats, savedCounts, vennPw);
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
						printCounts(feats, savedCounts, vennPw);
					}
				}
			}
		}
		
		vennPw.close();
	}

	protected int printCounts(List<String> featTypes, Map<String,Integer> savedCounts, PrintWriter vennPw) 
	throws Exception
	{
		StringBuffer sb = new StringBuffer(5000);
		int nFeats = featTypes.size();
		
		// And the querier params
		MethylDbQuerier params = new MethylDbQuerier();
		params.methylTablePrefix = this.methylDbPrefix;
		params.setMinCTreads(this.minCTreads);
		params.setUseNonconversionFilter(!this.noNonconvFilter);
		params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);

		// Feat filts
		params.addFeatFilters(featTypes);
		
		ListUtils.setDelim("+");
		String featStr = ListUtils.excelLine(featTypes);
		
		int count = 0;
		
		for (String chr : Arrays.asList("chr11","chr12")) //MethylDbUtils.CHROMS) //  
		{		
			
			// Iterator uses DB connection and can use a ton of memory because
			// it loads all rows at once.  This stuff should really be added to iterator
			// class, but until it is , just iterate here over the chromosome
			MAXCOORD = GoldAssembly.chromLengthStatic(chr, "hg18");
			for (int c = 0; c < MAXCOORD; c += STEP)
			{
				int last = Math.min(c+STEP-1, MAXCOORD);
				params.clearRangeFilters();
				params.addRangeFilter(chr, c, last);
//				params.addRangeFilter(chr);
			
				CpgIterator cpgit = new CpgIterator(params);
				int numRows = cpgit.getCurNumRows();
				count += numRows;
			}
		}
		

		if ((nFeats > 1) && (vennPw != null))
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
			double divBy = MatUtils.nanMax(allVals) / 100.0;

			
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

			vennPw.printf("<H3>%s</H3>\n", featStr);
			vennPw.printf("<IMG SRC=\"%s\">\n", url);
		}

        
		// Output
		sb.append(String.format("%d,%s",count,featStr));
		for (int i = 0; i < nFeats; i++)
		{
			Integer singleCount = savedCounts.get(featTypes.get(i));
			if (singleCount != null)
			{
				sb.append(String.format(",%.2f%%", 100.0*(double)count/(double)singleCount));
			}
		}
		
		sb.append("\n");
		savedCounts.put(featStr, count);
		System.err.println("Setting " + featStr + ": " + count);
		System.out.append(sb);
		return count;
	}
	


	
}

