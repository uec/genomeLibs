package edu.usc.epigenome.scripts;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
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
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairs;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsBinnedAutocorr;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsPearsonAutocorr;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;



public class MethylDbToAutocorr {

	public int STEP = (int)1E6;
	public int MINCOORD = 0;
	public int MAXCOORD = (int)2.8E8;
	
	
	
	private static final String C_USAGE = "Use: MethylDbToAutocorr -range 0.0 -range 0.5 -range 1.0 -outPrefix out -table methylCGsRich_normal010310_ -table methylCGsRich_tumor011010_ " + 
		" -windSize 2000 -minCTreads 2 -maxOppStrandAfrac 0.10 -maxNextNonGfrac 0.10 -noNonconvFilter";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-table",usage="Prefix for DB table (default methylCGsRich_normal010310_)")
    protected List<String> tables = new ArrayList<String>();
    @Option(name="-usePearson",usage="Does not stratify by bins.  Uses simple correlation coefficient metric")
    protected boolean usePearson = false;
    @Option(name="-sameStrand",usage="Only count pairs on the same strand (default false)")
    protected boolean sameStrand = false;
    @Option(name="-range",usage="All the range breakpoints that you want (default 0.0, 0.33, 0.67, 1.0)")
    protected List<Double> range = null;
    @Option(name="-withinFeat",usage="A featType from the features table")
    protected String withinFeat = null;
    @Option(name="-featFlank",usage="Flank size to use with the feature")
    protected int featFlank = 0;
    @Option(name="-outPrefix",usage="Output files will have this name")
    protected String outPrefix = "wiggleTester";
    @Option(name="-windSize",usage="starting window size (500)")
    protected int windSize = 2000;
    @Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
    protected int minCTreads = 2;
    @Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
    		"This sets a maximum number of observed As on the opposite strand (default 0.1)")
    protected double maxOppStrandAfrac = 0.1;
    @Option(name="-maxNextNonGfrac",usage="If the base following the C has more than this ratio of non-G bases, we don't count it. (default 0.1)")
    protected double maxNextNonGfrac = 0.1;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToAutocorr().doMain(args);
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

			if( tables.size() < 1)
			{
				System.err.println("Must specify at least one table");
				System.err.println(C_USAGE);
				System.exit(1);
			}

			if( arguments.size() != 0 ) {
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
		
		if (range == null)
		{
			range = Arrays.asList(0.0,0.3333,0.6667,1.0);
		}
		else
		{
			Collections.sort(range);
		}
		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.SEVERE);
		
		int nTables = this.tables.size();
		
//		MINCOORD = 7000000;
//		MAXCOORD = 10000000;
//		STEP = MAXCOORD-MINCOORD+1;;
		

		// Setup output files
		List<PrintWriter> pws = new ArrayList<PrintWriter>();
		for (int i = 0; i < nTables; i++)
		{
			// Start table
			String tab = this.tables.get(i);

			String strandSec = (this.sameStrand) ? ".sameStrand" : "";
			String featSec = (this.withinFeat==null) ? "" : String.format(".%s-f%d", this.withinFeat, this.featFlank);
			String outFn = String.format("Autocorr.%s.%s%s%s.wind%d.csv", 
					this.outPrefix, tab, featSec, strandSec, this.windSize);
			PrintWriter pw = new PrintWriter(new FileOutputStream(outFn));
			pws.add(pw);
		}
	
		// Setup streaming params
		MethylDbQuerier params = new MethylDbQuerier();
		params.setMinCTreads(this.minCTreads);
		params.setUseNonconversionFilter(!this.noNonconvFilter);
		params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
		params.setMaxNextNonGfrac(this.maxNextNonGfrac);
		if (this.withinFeat!=null) params.addFeatFilter(this.withinFeat, this.featFlank);

		
		// If we are doing correlation, we need the mean and SD.
		CpgMethLevelSummarizer[] methSumms = getMethSummarizers(params);
		
		
		// Setup autocorr objects
		CpgWalkerParams walkerParams = new CpgWalkerParams();
		walkerParams.maxWindSize = this.windSize;
		walkerParams.minWindSize = this.windSize;
		walkerParams.minCpgs = 2;
		CpgWalkerAllpairs[] autocorrs = new CpgWalkerAllpairs[nTables];
		for (int i = 0; i < nTables; i++)
		{
			// Autocorr counter
			if (this.usePearson)
			{
				double mean = methSumms[i].getValMean(false);
				double sd = methSumms[i].getValStdev();
				System.err.printf("mean=%.3f, sd=%.3f (table %s)\n",mean,sd,this.tables.get(i));
				autocorrs[i] = new CpgWalkerAllpairsPearsonAutocorr(walkerParams, sameStrand, mean, sd);
			}
			else
			{
				autocorrs[i] = new CpgWalkerAllpairsBinnedAutocorr(walkerParams, sameStrand, this.range);
			}
	
			// Header
			PrintWriter pw = pws.get(i);
			String header = autocorrs[i].headerStr();
			if (header != null) pw.println();
		}
		



		for (String chr :  Arrays.asList("chr11"))//MethylDbUtils.CHROMS) // 
		{
		

			// Iterator uses DB connection and can use a ton of memory because
			// it loads all rows at once.  This stuff should really be added to iterator
			// class, but until it is , just iterate here over the chromosome
			int onCpg = 0;

			
			for (int c = MINCOORD; c < MAXCOORD; c += STEP)
			{
				System.err.printf("LOADING NEW WIND: %d-%d\n",c,c+STEP-1);

				params.clearRangeFilters();
				params.addRangeFilter(chr, c, c+STEP-1); // *** TESTING , REPLACE ****
				CpgIteratorMultisample cpgit = new CpgIteratorMultisample(params, this.tables);
				//int numCpgs = cpgit.getCurNumRows();
				while (cpgit.hasNext())
				{
					Cpg[] cpgs = cpgit.next();

					// Stream Cpgs
					for (int i = 0; i < nTables; i++)
					{
						autocorrs[i].streamCpg(cpgs[i]);
					}

					//
					if ((onCpg % 1E5)==0) System.err.printf("On Cpg #%d, domain %s\n", onCpg, autocorrs[0].windStr());
					onCpg++;
				}


			}
		}
	
		// Print autocorrs and close
		for (int i = 0; i < nTables; i++)
		{
			PrintWriter pw = pws.get(i);
			CpgWalkerAllpairs autocorr = autocorrs[i];
			
			pw.append(autocorr.toCsvStr());
			
			pw.close();
		}
		
		
	} // Main

	private CpgMethLevelSummarizer[] getMethSummarizers(MethylDbQuerier params) 
	throws Exception
	{

		System.err.println("Getting mean/sd");
		
		int nTabs = this.tables.size();
		CpgMethLevelSummarizer[] summs = new CpgMethLevelSummarizer[nTabs];
		for (int i=0; i<nTabs; i++)
		{
			summs[i] = new CpgMethLevelSummarizer(params);
		}
		
		for (String chr : Arrays.asList("chr11")) //  MethylDbUtils.CHROMS) //
		{
			// Iterator uses DB connection and can use a ton of memory because
			// it loads all rows at once.  This stuff should really be added to iterator
			// class, but until it is , just iterate here over the chromosome
			int onCpg = 0;

			for (int c = MINCOORD; c < MAXCOORD; c += STEP)
			{
				System.err.printf("LOADING NEW WIND (MEAN/SD pass): %d-%d\n",c,c+STEP-1);

				params.clearRangeFilters();
				params.addRangeFilter(chr, c, c+STEP-1); 
				CpgIteratorMultisample cpgit = new CpgIteratorMultisample(params, this.tables);
				while (cpgit.hasNext())
				{
					Cpg[] cpgs = cpgit.next();

					// Stream Cpgs
					for (int i = 0; i < nTabs; i++)
					{
						summs[i].streamCpg(cpgs[i]);
					}

					//
					onCpg++;
				}
			}
		}
		
		return summs;
	}

	public static double safeVal(double v)
	{
		return ((Double.isNaN(v) || Double.isInfinite(v)) ? 0.0 : v);
	}
	
}


