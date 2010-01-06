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



public class MethylDbToSummaryStats {

	private static final String C_USAGE = "Use: MethylDbToSummaryStats table1prefix table2prefix .. ";
	
	private static List<Color> colors = Arrays.asList(Color.BLUE, Color.RED, Color.GREEN, Color.ORANGE, Color.BEIGE, Color.BLACK, Color.AQUAMARINE);
	
//    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
//    protected boolean noNonconvFilter = false;
	// receives other command line parameters than options
	@Argument
	private List<String> tablePrefixes = new ArrayList<String>();
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToSummaryStats().doMain(args);
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
		
		// Start output
		StringWriter sw = new StringWriter(200000);
		PrintWriter pw = new PrintWriter(sw);
		

		// Go though features , then samples
		List<String> featTabs = MethylDbUtils.SUMMARY_FEATURES1;//FeatIterator.AllFeatTablePrefixes(); 
		int nS = tablePrefixes.size();
		int nF = featTabs.size();
		
		// Go through features than samples
		
		System.err.printf("nS=%d, nF=%d\n", nS,nF);
		CpgSummarizer[][][] summarizerLists = new CpgSummarizer[nF][][];
		int numRows = 0;
		for (int i = 0; i < nF; i++)
		{
			String feat = featTabs.get(i);
			summarizerLists[i] = new CpgSummarizer[nS][];
			for (int j = 0; j < nS; j++)
			{
				String sample = tablePrefixes.get(j);
				System.err.printf("Filling summarizerRow[%d]\n", numRows);
				summarizerLists[i][j] = summarizers(sample, feat);
			}
		}
		
		// Make charts
		int nSumms = summarizerLists[0][0].length;
		for (int s = 0; s<nSumms; s++)
		{
			double natMin = summarizerLists[0][0][s].getNaturalMin();
			double natMax = summarizerLists[0][0][s].getNaturalMax();
			
			// Each of the samples should be data class
			List<Plot> samplePlots = new ArrayList<Plot>(nS);
			for (int j = 0; j < nS; j++)
			{
				// Get the data values
				double[] valMeans = new double[nF];
				for (int i =0; i<nF; i++)
				{
					valMeans[i] = summarizerLists[i][j][s].getValMean();
				}
				System.err.printf("Data series for %s has %d samples.  natmin=%f, natmax=%f\n",
						summarizerLists[0][j][s].getSampleName(), valMeans.length, natMin, natMax);
				MatUtils.nansToVal(valMeans, -0.1);
				Data data= DataUtil.scaleWithinRange(natMin, natMax, valMeans);
				//Data data= DataUtil.scale(valMeans);
				Color c = colors.get(j);
				
				String sampleName = summarizerLists[0][j][s].getSampleName();
				sampleName = sampleName.replaceAll("methylCGsRich_", "");
				BarChartPlot plot = Plots.newBarChartPlot(data, c, sampleName);
				samplePlots.add(plot);
			}
			
			// Chart type
			String summarizerLabel = summarizerLists[0][0][s].getLabel();
			
			BarChart chart = GCharts.newBarChart(samplePlots);
			chart.setSize(950, 200); // Make sure to make it wide enough, it doesn't autoscale
	        chart.setDataEncoding(DataEncoding.SIMPLE);

	        AxisLabels xAxis = AxisLabelsFactory.newNumericRangeAxisLabels(1, nF);
	       	chart.addXAxisLabels(xAxis);
//	       	chart.addXAxisLabels(AxisLabelsFactory.newAxisLabels("Features"));
	       	
	        AxisLabels yAxis = AxisLabelsFactory.newNumericRangeAxisLabels(natMin, natMax);
	       	chart.addYAxisLabels(yAxis);
	       	chart.addRightAxisLabels(yAxis);
	
	       	chart.setTitle(summarizerLabel);
	       	
	       	
	        
	        pw.append("<P>");
			pw.append("<IMG ");
			pw.append(" SRC=\"");
			pw.append(chart.toURLString());
			pw.append("\" ALT=\"google chart\">\n");
			pw.append("</P>\n");
		}
		
		pw.println(CpgSummarizer.htmlTableStart());
		
		for (int s = 0; s<nSumms; s++)
		{
			for (int i =0; i<nF; i++)
			{
				for (int j = 0; j < nS; j++)
				{
					//System.err.printf("summarizerRows[%d][%d]\n", i,j);
					pw.println(summarizerLists[i][j][s].htmlTableRow());
				}
			}
		}
	
		
		pw.println(CpgSummarizer.htmlTableFinish());
		
		
		System.out.print(sw.toString());
	}
	
	
	public CpgSummarizer[] summarizers(String sampleTablePrefix, String feat)
	throws Exception
	{

		// Don't filter
		MethylDbQuerier querier = new MethylDbQuerier();
		querier.setMinCTreads(0);
		querier.setUseNonconversionFilter(false);
		querier.setMaxOppstrandAfrac(Double.POSITIVE_INFINITY);

		querier.setMethylTablePrefix(sampleTablePrefix);
		if (feat!=null) querier.addFeatFilter(feat);

		// Setup summarizers
		CpgDensitySummarizer densSummarizer = new CpgDensitySummarizer(querier);
		List<CpgSummarizer> summs = new ArrayList<CpgSummarizer>(10);
		summs.add(new CpgMethLevelSummarizer(querier));
		summs.add(densSummarizer);
		summs.add(new CpgNonconversionSummarizer(querier));
		summs.add(new CpgCoverageSummarizer(querier));
		summs.add(new CpgDeaminationSummarizer(querier));

		for (String chrom : MethylDbUtils.CHROMS)
		{
			querier.clearRangeFilters();
			querier.addRangeFilter(chrom);

//			querier.clearRangeFilters();
//			querier.addRangeFilter(chrom,17322212,17877427);
			
			densSummarizer.addRanges(querier);
			
			CpgIterator cpgit = new CpgIterator(querier);
			int numRows = cpgit.getCurNumRows();
			Logger.getAnonymousLogger().severe("Found " + numRows + " cpgs");
			int i = 0;
			while (cpgit.hasNext())
			{
				Cpg cpg = cpgit.next();
				i++;
				if ((i%1E5)==0) System.err.println("On " + i + "/" + numRows + ":\t" + cpg.toString());

				// Stream cpgs
				for (CpgSummarizer summ : summs)
				{
					summ.streamCpg(cpg);
				}
			}
		}

//		String[] out = new String[summs.size()];
//		for (int i=0; i < summs.size(); i++)
//		{
//			CpgSummarizer summ = summs.get(i);
//			out[i] = summ.htmlTableRow();
//		}
		
		return summs.toArray(new CpgSummarizer[1]);
	}
	
}

