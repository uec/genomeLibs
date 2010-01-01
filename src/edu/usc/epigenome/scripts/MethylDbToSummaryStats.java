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

import com.sun.org.apache.xalan.internal.xsltc.compiler.Pattern;
import com.sun.tools.javac.code.Attribute.Array;

import edu.usc.epigenome.genomeLibs.FeatDb.FeatIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgCoverageSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgMethLevelSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;



public class MethylDbToSummaryStats {

	private static final String C_USAGE = "Use: MethylDbToSummaryStats table1prefix table2prefix .. ";
	
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
		pw.println(CpgSummarizer.htmlTableStart());
		

		// Go though features , then samples
		List<String> featTabs = FeatIterator.AllFeatTablePrefixes(); 
		int nS = tablePrefixes.size();
		int nF = featTabs.size();
		
		// Go through features than samples
		
		System.err.printf("nS=%d, nF=%d\n", nS,nF);
		String[][] summarizerRows = new String[nS*nF][]; // Each one will be a list of summarizers
		int numRows = 0;
		for (int i = 0; i < nF; i++)
		{
			String feat = featTabs.get(i);
			for (int j = 0; j < nS; j++)
			{
				String sample = tablePrefixes.get(j);
				System.err.printf("Filling summarizerRow[%d]\n", numRows);
				summarizerRows[numRows++] = rows(sample, feat);
			}
		}
		
		for (int j = 0; j < summarizerRows[0].length; j++)
		{
			for (int i = 0; i < numRows; i++)
			{
				System.err.printf("summarizerRows[%d][%d]\n", i,j);
				pw.println(summarizerRows[i][j]);
			}
		}
		
		pw.println(CpgSummarizer.htmlTableFinish());
		
		
		System.out.print(sw.toString());
	}
	
	
	public String[] rows(String sampleTablePrefix, String feat)
	throws Exception
	{

		// Don't filter
		MethylDbQuerier querier = new MethylDbQuerier();
		querier.setMinCTreads(0);
		querier.setUseNonconversionFilter(true);
		querier.setMaxOppstrandAfrac(0);

		querier.setMethylTablePrefix(sampleTablePrefix);
		if (feat!=null) querier.addFeatFilter(feat);

		// Setup summarizers
		List<CpgSummarizer> summs = new ArrayList<CpgSummarizer>(10);
		summs.add(new CpgMethLevelSummarizer(querier));
		summs.add(new CpgCoverageSummarizer(querier));

		for (String chrom : MethylDbUtils.CHROMS11)
		{
			querier.clearRangeFilters();
//			querier.addRangeFilter(chrom,17322212,17877427);
			querier.addRangeFilter(chrom);

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

		String[] out = new String[summs.size()];
		for (int i=0; i < summs.size(); i++)
		{
			CpgSummarizer summ = summs.get(i);
			out[i] = summ.htmlTableRow();
		}
		
		return out;
	}
	
}

