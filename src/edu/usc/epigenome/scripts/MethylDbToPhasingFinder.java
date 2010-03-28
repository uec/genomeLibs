package edu.usc.epigenome.scripts;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
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
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerPhasingFinder;



public class MethylDbToPhasingFinder {

	public int STEP = (int)1E6;
	public int MINCOORD = 0;
	public int MAXCOORD = (int)2.8E8;
	
	
	
	private static final String C_USAGE = "Use: MethylDbToPhasingFinder -outPrefix out -table methylCGsRich_normal010310_ -table methylCGsRich_tumor011010_ " + 
		" -windSize 2000 -minCTreads 2 -maxOppStrandAfrac 0.10 -maxNextNonGfrac 0.10 -noNonconvFilter";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-table",usage="Prefix for DB table (default methylCGsRich_normal010310_)")
    protected List<String> tables = new ArrayList<String>();
    @Option(name="-sameStrand",usage="Only count pairs on the same strand (default false)")
    protected boolean sameStrand = false;
    @Option(name="-matlabStyle",usage="MATLAB instead of wig output (default false)")
    protected boolean matlabStyle = false;
    @Option(name="-withinFeat",usage="A featType from the features table")
    protected List<String> withinFeats = new ArrayList<String>(10);
    @Option(name="-featFlank",usage="Flank size to use with the feature")
    protected int featFlank = 0;
    @Option(name="-nucFlank",usage="Flank size from nucleosome position (default 20)")
    protected int nucFlank = 20;
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
		new MethylDbToPhasingFinder().doMain(args);
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
		
		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.SEVERE);
		
		int nTables = this.tables.size();
		

		
		
		// Loop through within feats.  If more than one is specified, throw in a null one for fun.  Also add
		// it if we haven't specified any
		if (withinFeats.size() != 1) withinFeats.add(null);

		for (String withinFeat : withinFeats)
		{

			// Setup output files
			List<PrintWriter> pws = new ArrayList<PrintWriter>();
			PrintWriter pwArr[][] = new PrintWriter[nTables][2];
			for (int i = 0; i < nTables; i++)
			{
				// Start table
				String tab = this.tables.get(i);

				for (int j = 0; j <= 1; j++)
				{
					String extension = (this.matlabStyle) ? ".csv" : ".wig";
					String strandSec = (this.sameStrand) ? ".sameStrand" : "";
					String featSec = String.format(".%s-f%d", (withinFeat==null)?"all":withinFeat, this.featFlank);
					String outFn = String.format("Phasing.%s.%s%s%s.wind%d.nucflank%d.%d%s", 
							this.outPrefix, tab, featSec, strandSec, this.windSize, this.nucFlank, j,extension);
					PrintWriter pw = new PrintWriter(new FileOutputStream(outFn));
					String name = String.format("%s %s", (j==0)?"half Nuc":"full Nuc", tab);
					if (!this.matlabStyle) pw.printf("track type=wiggle_0 name=\"%s\" description=\"%s\"\n",name,name);	
					pwArr[i][j] = pw;
					pws.add(pw);
				}
			}

			// Setup streaming params
			MethylDbQuerier params = new MethylDbQuerier();
			params.setMinCTreads(this.minCTreads);
			params.setUseNonconversionFilter(!this.noNonconvFilter);
			params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
			params.setMaxNextNonGfrac(this.maxNextNonGfrac);


			// Autocorr params needed for pass 1 and 2
			CpgWalkerParams walkerParams = new CpgWalkerParams();
			walkerParams.maxScanningWindSize = this.windSize;
			walkerParams.minScanningWindSize = this.windSize;
			walkerParams.minScanningWindCpgs = 2;

			
			// For the background, it makes a big difference whether we add this filter
			// before or after doing the 1st pass to calculate mean/sd.
			if (withinFeat!=null) params.addFeatFilter(withinFeat, this.featFlank);



			CpgWalkerPhasingFinder[] finders = new CpgWalkerPhasingFinder[nTables];
			for (int i = 0; i < nTables; i++)
			{
				//PrintWriter pw = pws.get(i);
				finders[i] = new CpgWalkerPhasingFinder(walkerParams, sameStrand, 
						185, pwArr[i][0],pwArr[i][1], this.matlabStyle, this.nucFlank);
			}

			// And run!
			this.streamFinders(params, walkerParams, finders);



			// Print autocorrs and close
			for (PrintWriter pw : pws)
			{
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
	private void streamFinders(MethylDbQuerier params, CpgWalkerParams walkerParams,
			CpgWalkerPhasingFinder[] finders) 
	throws Exception
	{
		int nTabs = this.tables.size();
			

		for (String chr : MethylDbUtils.CHROMS) //  Arrays.asList("chr11","chr12","chr13")) //
		{

			// Stream Cpgs
			for (int i = 0; i < nTabs; i++)
			{
				finders[i].setCurChr(chr);
			}

			// Iterator uses DB connection and can use a ton of memory because
			// it loads all rows at once.  This stuff should really be added to iterator
			// class, but until it is , just iterate here over the chromosome
			int onCpg = 0;

//			MINCOORD = 7000000;
//			MAXCOORD = 8000000;
//			STEP = Math.min(STEP, MAXCOORD-MINCOORD+1);
			

			for (int c = MINCOORD; c < MAXCOORD; c += STEP)
			{
				System.err.printf("LOADING NEW WIND: %d-%d\n",c,c+STEP-1);

				params.clearRangeFilters();
				params.addRangeFilter(chr, c, c+STEP-1);
				CpgIteratorMultisample cpgit = new CpgIteratorMultisample(params, this.tables);
				//int numCpgs = cpgit.getCurNumRows();
				while (cpgit.hasNext())
				{
					Cpg[] cpgs = cpgit.next();

					// Stream Cpgs
					for (int i = 0; i < nTabs; i++)
					{
						//System.err.printf("Streaming cpg at pos %d to table %s\n",cpgs[i].chromPos,this.tables.get(i));
						finders[i].streamCpg(cpgs[i]);
					}

					//
					if ((onCpg % 1E5)==0) System.err.printf("On Cpg #%d, domain %s\n", onCpg, finders[0].windStr());
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


