package edu.usc.epigenome.scripts;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;


import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


public class MethylDbToMultiMeth {

	private static final String C_USAGE = "Use: MethylDbToMultiMeth -chrom chr1 -chrom chr 2 tablePrefix tablePrefix2 .. > stats.csv";
	
    @Option(name="-minCt",usage="minimum number of C or T reads (in all samples) (default 1)")
    protected int minCt = 1;
    @Option(name="-chrom",usage="A chromosome to use")
    protected List<String> chroms = null;
	// receives other command line parameters than options
	@Argument
	private List<String> tablePrefixes = new ArrayList<String>();
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToMultiMeth().doMain(args);
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
		

		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.FINE);
		
		int nS = tablePrefixes.size();
		
		// Go through features than samples
		// Don't filter
		MethylDbQuerier querier = new MethylDbQuerier();
		querier.setMinCTreads(this.minCt);
		querier.setUseNonconversionFilter(true);
		querier.setMaxOppstrandAfrac(0.2);

		// Setup counters
		if (this.chroms == null) chroms = MethylDbUtils.CHROMS;
		for (String chrom : chroms) 
		{

			// Stupid JDBC tries to load entire chromosome into memory at once,
			// which is too much.
			final int STEP = (int)1E7;
			final int MAXCOORD = (int)2.8E8;
			for (int c = 0; c < MAXCOORD; c += STEP)
			{
			
				Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Starting chrom " + chrom + " from " + (c+1) + " to " + (c+STEP) + "\n");
		
			querier.clearRangeFilters();
			querier.addRangeFilter(chrom, c+1, c+STEP);
			
			// This will use too much memory if we don't set "allowNumRows" to false
			CpgIteratorMultisample cpgit = new CpgIteratorMultisample(querier, tablePrefixes);
			int numRows = cpgit.getCurNumRows();
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Found " + numRows + " cpgs");
			int i = 0;
			while (cpgit.hasNext())
			{
				Cpg[] cpg = cpgit.next();
				i++;
				if ((i%1E5)==0) 
				{
					System.err.println("On " + i + "/" + numRows + ":\t" + cpg.toString());
				}

				for (int j = 0; j < nS; j++)
				{
					if (j>0) System.out.append(',');
					System.out.append(Double.toString(cpg[j].fracMeth(true)));
				}
				System.out.println();
				
			}
		}

		}

	}
	
	
}

