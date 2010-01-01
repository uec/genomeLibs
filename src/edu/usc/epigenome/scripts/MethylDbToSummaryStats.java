package edu.usc.epigenome.scripts;

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

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;



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
		
		
		MethylDbQuerier querier = new MethylDbQuerier();
		querier.setMinCTreads(4);
		querier.setUseNonconversionFilter(true);
		querier.setMaxOppstrandAfrac(0.2);
		querier.setMethylTablePrefix(tablePrefixes.get(0));
		
//		querier.addFeatFilter("hg18.ES.H3K27me3.HMM.gtf");
		querier.addRangeFilter("chr11", 30683099, 31960706);
		querier.addRangeFilter("chr11", 35683099, 36960706);
		
//		CpgIterator cpgit = new CpgIterator(querier);
//		int numRows = cpgit.getCurNumRows();
//		Logger.getAnonymousLogger().severe("Found " + numRows + " cpgs");
//		int i = 0;
//		while (cpgit.hasNext())
//		{
//			Cpg cpg = cpgit.next();
//			i++;
//			// System.err.println("On " + i + "/" + numRows + ":\t" + cpg.toString());
//		}
		
		CpgIteratorMultisample cpgit = new CpgIteratorMultisample(querier, Arrays.asList("methylCGsRich_normal123009_", "methylCGsRich_tumor123009_"));
		int numRows = cpgit.getCurNumRows();
		Logger.getAnonymousLogger().severe("Found " + numRows + " cpgs");
		int i = 0;
		while (cpgit.hasNext())
		{
			Cpg[] cpg = cpgit.next();
			i++;
			//System.err.println("On " + i + "/" + numRows + ":\t" + cpg[0].toString() + "\t" + cpg[1].toString());
		}

	}
	
}

