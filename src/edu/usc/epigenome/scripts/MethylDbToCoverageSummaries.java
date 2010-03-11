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


public class MethylDbToCoverageSummaries {

	private static final String C_USAGE = "Use: MethylDbToSummaryStats tablePrefix -withinFeat featA -withinFeat featB > stats.csv";
	
//    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
//    protected boolean noNonconvFilter = false;
	// receives other command line parameters than options
	@Option(name="-doubleStranded",usage="If set, counts include both strands (default false)")
	protected boolean doubleStranded = false;
	@Option(name="-withinFeat",usage="make separate count listings ones for each feature type")
	protected List<String> withinFeats = new ArrayList<String>(20);
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
		

		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.FINE);
		
		int nS = tablePrefixes.size();
		
		// Go through features than samples
		// Don't filter
		MethylDbQuerier querier = new MethylDbQuerier();
		querier.setMinCTreads(0);
		querier.setMaxNextNonGfrac(1.0);
		querier.setUseNonconversionFilter(false);
		querier.setMaxOppstrandAfrac(Double.POSITIVE_INFINITY);

		querier.setMethylTablePrefix(tablePrefixes.get(0));

		// One for each feature.  Add one for global
		withinFeats.add("global");
		int nF = withinFeats.size();
		
		
		// Setup counters
		long totalUniqueCpgs = 0;
		long totalMeasurements = 0;
		SortedMap<Integer,Integer>[] featCounts = new SortedMap[nF];
		for (int i = 0; i < nF; i++)
		{
			featCounts[i] = new TreeMap<Integer,Integer>();
		}
			
		for (String chrom : MethylDbUtils.CHROMS) // MethylDbUtils.TEST_CHROMS) // 
		{
			// Stupid JDBC tries to load entire chromosome into memory at once,
			// which is too much.
			int STEP = (int)1E6;
			int MINCOORD = 0;
			int MAXCOORD = (int)2.8E8;
			
//			MAXCOORD = 10000000;
//			STEP = Math.min(STEP, MAXCOORD-MINCOORD+1);

			for (int c = MINCOORD; c < MAXCOORD; c += STEP)
			{

				Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Starting chrom " + chrom + " from " + c+1 + " to " + c+STEP + "\n");

				querier.clearRangeFilters();
				querier.addRangeFilter(chrom, c+1, c+STEP);

				for (int f = 0; f < nF; f++)
				{

					String withinFeat = withinFeats.get(f);
					boolean global = withinFeat.matches("global");
					SortedMap<Integer,Integer> counts = featCounts[f];

					// Add feature filter
					querier.clearFeatFilters();
					if (!global) querier.addFeatFilter(withinFeat);


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
							System.err.printf("\tTree has %d items\n", counts.size());
						}

						int count = cpg[0].totalReads;
						if (this.doubleStranded) count += cpg[0].totalReadsOpposite;

						if (global)
						{
							if (count >= 1) totalUniqueCpgs++;
							totalMeasurements += count;
						}


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
			}
		}

		// Output
		//System.out.printf("Total unique CpG=%d, total CpG measurements=%d\n",totalUniqueCpgs, totalMeasurements);
		int maxCvg = 0;
		for (int f = 0; f < nF; f++)
		{
			maxCvg = Math.max(maxCvg, featCounts[f].lastKey());
		}
		
		// Header
		System.out.print("coverage");
		for (int f = 0; f < nF; f++)
		{
			System.out.printf(",%s",withinFeats.get(f));
		}
		System.out.println();
		
		
		for (int i = 0; i <= maxCvg; i++)
		{
			System.out.printf("%d",i);
			for (int f = 0; f < nF; f++)
			{
				Integer countObj = featCounts[f].get(new Integer(i));
				int count = (countObj==null) ? 0 : countObj.intValue();
				System.out.printf(",%d",count);
			}
			System.out.println();
		}
		

	}
	
	
}

