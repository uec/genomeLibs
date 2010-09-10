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

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairs;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsAutocorrByread;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsBinnedAutocorr;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsPearsonAutocorr;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;



public class MethylDbToAutocorrByread {

	public int STEP = (int)1E6;
	public int MINCOORD = 0;
	public int MAXCOORD = (int)2.8E8;
	
	
	
	private static final String C_USAGE = "Use: MethylDbToAutocorrByread -outPrefix out -table methylCGsRich_normal010310_ " + 
		" -windSize 2000 -minCTreads 2 -maxOppStrandAfrac 0.10 -maxNextNonGfrac 0.10 -noNonconvFilter";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-table",usage="Prefix for DB table (default methylCGsRich_normal010310_)")
    protected String table = null;
//    @Option(name="-usePearson",usage="Does not stratify by bins.  Uses simple correlation coefficient metric")
//    protected boolean usePearson = false;
//    @Option(name="-sameStrand",usage="Only count pairs on the same strand (default false)")
//    protected boolean sameStrand = false;
//    @Option(name="-range",usage="All the range breakpoints that you want (default 0.0, 0.33, 0.67, 1.0)")
//    protected List<Double> range = null;
    @Option(name="-withinFeat",usage="A featType from the features table")
    protected List<String> withinFeats = new ArrayList<String>(10);
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
		new MethylDbToAutocorrByread().doMain(args);
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

			if( table == null)
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
		

		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
		
		int nTables = 1;
		

		
		
		// Loop through within feats.  If more than one is specified, throw in a null one for fun.  Also add
		// it if we haven't specified any
		if (withinFeats.size() != 1) withinFeats.add(null);

		for (String withinFeat : withinFeats)
		{



			// Setup streaming params
			MethylDbQuerier params = new MethylDbQuerier();
			params.setMethylTablePrefix(this.table);
			params.setReadTableSchema();
			params.setMinCTreads(this.minCTreads);
			params.setUseNonconversionFilter(!this.noNonconvFilter);
			params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
			params.setMaxNextNonGfrac(this.maxNextNonGfrac);


			// Autocorr params needed for pass 1 and 2
			CpgWalkerParams walkerParams = new CpgWalkerParams();
			walkerParams.maxScanningWindSize = this.windSize;
			walkerParams.minScanningWindSize = this.windSize;
			walkerParams.minScanningWindCpgs = 2;
			walkerParams.methylParams = params;

			
			// For the background, it makes a big difference whether we add this filter
			// before or after doing the 1st pass to calculate mean/sd.
			if (withinFeat!=null) params.addFeatFilter(withinFeat, this.featFlank);


			// Setup output files and autocorrs
			List<CpgWalkerAllpairs> autocorrs = new ArrayList<CpgWalkerAllpairs>(10);
			List<PrintWriter> pws = new ArrayList<PrintWriter>(10);
			
			String tab = this.table;
			String featSec = String.format(".%s-f%d", (withinFeat==null)?"all":withinFeat, this.featFlank);
			String outFnPrefix = String.format("Autocorr.%s.%s%s.wind%d", 
					this.outPrefix, tab, featSec, this.windSize);
			
			final int N_CONDITIONS = 5;
			final int START_COND = 1;
			final int END_COND = 5;

//			final int N_CONDITIONS = 1;
//			final int START_COND = 3;
//			final int END_COND = 3;
			
			for (int i = START_COND; i <= END_COND; i++)
			{
				String typeStr = null;
				boolean sameStrand = false;
				boolean oppStrand = false;
				boolean sameRead = false;
				boolean differentRead = false;
				
				switch(i)
				{
				case 1: typeStr = "anyStrandAnyRead" ; sameStrand=false; oppStrand=false; sameRead=false; differentRead=false; break;
				case 2: typeStr = "oppStrandAnyRead" ; sameStrand=false; oppStrand=true; sameRead=false; differentRead=false; break;
				case 3: typeStr = "sameStrandAnyRead" ; sameStrand=true; oppStrand=false; sameRead=false; differentRead=false; break;
				case 4: typeStr = "sameStrandSameRead" ; sameStrand=true; oppStrand=false; sameRead=true; differentRead=false; break;
				case 5: typeStr = "sameStrandDifferentRead" ; sameStrand=true; oppStrand=false; sameRead=false; differentRead=true; break;
				}
				
				PrintWriter pw = new PrintWriter(new FileOutputStream(outFnPrefix + "-" + typeStr + ".csv"));
				pws.add(pw);
				CpgWalkerAllpairs autocorr = new CpgWalkerAllpairsAutocorrByread(walkerParams, sameStrand, oppStrand, sameRead, differentRead);
				autocorrs.add(autocorr);
				pw.println(autocorr.headerStr());
			}
			
			// And run
			this.streamAutocorrs(params, walkerParams, autocorrs);



			// Print autocorrs and close
			for (int i = 0; i < N_CONDITIONS; i++)
			{
				PrintWriter pw = pws.get(i);
				CpgWalkerAllpairs autocorr = autocorrs.get(i);
				pw.append(autocorr.toCsvStr());
				pw.close();
			}

		}

		
	} // Main

	
	/**
	 * @param params
	 * @param walkerParams
	 * @param autocorrs If passed in , we use these ones.  If null, we create Pearson autocorrs from scratch,
	 *        and run them using mean 0 SD 0 (useful for doing a first pass to calculate mean/sd for second pass)
	 * @return
	 * @throws Exception
	 */
	private void streamAutocorrs(MethylDbQuerier params, CpgWalkerParams walkerParams,
			List<CpgWalkerAllpairs> autocorrs) 
	throws Exception
	{


		for (String chr : Arrays.asList("chr2")) // MethylDbUtils.CHROMS) //  
		{


			// Iterator uses DB connection and can use a ton of memory because
			// it loads all rows at once.  This stuff should really be added to iterator
			// class, but until it is , just iterate here over the chromosome
			int onCpg = 0;

//			MINCOORD = 7000000;
//			MAXCOORD = 7010000;
//			STEP = Math.min(STEP, MAXCOORD-MINCOORD+1);
			

			for (int c = MINCOORD; c < MAXCOORD; c += STEP)
			{
				System.err.printf("LOADING NEW WIND: %d-%d\n",c,c+STEP-1);

				params.clearRangeFilters();
				params.addRangeFilter(chr, c, c+STEP-1);
				CpgIterator cpgit = new CpgIterator(params);
				//int numCpgs = cpgit.getCurNumRows();
				while (cpgit.hasNext())
				{
					Cpg cpg = cpgit.next();
					
					for (CpgWalkerAllpairs autocorr : autocorrs)
					{
						autocorr.streamCpg(cpg);
					}

					CpgWalkerAllpairs firstAutocorr = autocorrs.get(0);
					
					if ((onCpg % 1E5)==0) System.err.printf("On Cpg #%d, domain %s\n", onCpg, firstAutocorr.windStr());
					onCpg++;
				}


			}
		}
	}
	

	public static double safeVal(double v)
	{
		return ((Double.isNaN(v) || Double.isInfinite(v)) ? 0.0 : v);
	}
	
}


