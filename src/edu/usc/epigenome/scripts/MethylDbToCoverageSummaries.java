package edu.usc.epigenome.scripts;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
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
import com.googlecode.charts4j.GCharts;
import com.googlecode.charts4j.Plot;
import com.googlecode.charts4j.Plots;
import com.sun.org.apache.xalan.internal.xsltc.compiler.Pattern;
import com.sun.tools.javac.code.Attribute.Array;

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



public class MethylDbToCoverageSummaries {

	private static final String C_USAGE = "Use: MethylDbToSummaryStats tablePrefix > stats.csv";
	
	private static List<Color> colors = Arrays.asList(Color.BLUE, Color.RED, Color.GREEN, Color.ORANGE, Color.BEIGE, Color.BLACK, Color.AQUAMARINE);
	
//    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
//    protected boolean noNonconvFilter = false;
	// receives other command line parameters than options
	@Argument
	private List<String> tablePrefixes = new ArrayList<String>();
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToCoverageSummaries().doMain(args);
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

			if(tablePrefixes.size() < 1 ) {
				System.err.println(C_USAGE);
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
		

		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.SEVERE);
		
		int nS = tablePrefixes.size();
		
		// Go through features than samples
		// Don't filter
		MethylDbQuerier querier = new MethylDbQuerier();
		querier.setMinCTreads(1);
		querier.setUseNonconversionFilter(false);
		querier.setMaxOppstrandAfrac(Double.POSITIVE_INFINITY);

		querier.setMethylTablePrefix(tablePrefixes.get(0));

		// Setup counters
		long totalUniqueCpgs = 0;
		long totalMeasurements = 0;
		SortedMap<Integer,Integer> counts = new TreeMap<Integer,Integer>();

		for (String chrom :Arrays.asList("chr11"))
		{
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Starting chrom " + chrom + "\n");
		
			querier.clearRangeFilters();
			querier.addRangeFilter(chrom, 1000000, 1050000);
			
			CpgIterator cpgit = new CpgIterator(querier);
			int numRows = cpgit.getCurNumRows();
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Found " + numRows + " cpgs");
			int i = 0;
			while (cpgit.hasNext())
			{
				Cpg cpg = cpgit.next();
				i++;
				if ((i%1E5)==0) System.err.println("On " + i + "/" + numRows + ":\t" + cpg.toString());

				int count = cpg.totalReads + cpg.totalReadsOpposite;
				totalUniqueCpgs++;
				totalMeasurements += count;
				
				
				Integer cvg = counts.get(count);
				if (cvg == null)
				{
					cvg = new Integer(0);
					counts.put(new Integer(count), cvg);
				}
				cvg = cvg + 1;
				//System.err.printf("New val for counts->{%d} = %d\n", count, cvg);
				counts.put(new Integer(count), cvg);
			}
		}
		
		// Output
		System.out.printf("Total unique CpG=%d, total CpG measurements=%d\n",totalUniqueCpgs, totalMeasurements);
		int maxCvg = counts.lastKey();
		for (int i = 0; i < maxCvg; i++)
		{
			Integer countObj = counts.get(new Integer(i));
			int count = (countObj==null) ? 0 : countObj.intValue();
			System.out.printf("%d,%d\n", i,count);
		}
		

	}
	
	
}

