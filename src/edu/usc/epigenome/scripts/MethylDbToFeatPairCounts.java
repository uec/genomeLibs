package edu.usc.epigenome.scripts;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
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
	
	
	// receives other command line parameters than options
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
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
		
		// Go through individual feats.
		for (String featType : featTypes)
		{
			List<String> feats = new ArrayList<String>(1);
			feats.add(featType);
			printCounts(feats);
		}
		
		// Then feat pairs
		for (int i = 0; i < nFeats; i++)
		{
			for (int j = i+1; j < nFeats; j++)
			{
				List<String> feats = new ArrayList<String>(2);
				feats.add(featTypes.get(i));
				feats.add(featTypes.get(j));
				printCounts(feats);
			}
		}
		

	}

	protected void printCounts(List<String> featTypes) 
	throws Exception
	{

		StringBuffer sb = new StringBuffer(5000);
		
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
		for (String chr : Arrays.asList("chr11")) //MethylDbUtils.CHROMS) //  
		{		
			params.clearRangeFilters();
			params.addRangeFilter(chr);
			
			CpgIterator cpgit = new CpgIterator(params);
			int numRows = cpgit.getCurNumRows();
			count += numRows;
			
			System.err.printf("%d, %s %s\n",numRows,chr,featStr);
		}
		
	}
	


	
}

