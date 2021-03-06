package edu.usc.epigenome.scripts;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.seq.StrandedFeature;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerDomainFinder;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerDomainFinderMethRange;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerMultisampleMethRange_DEPR;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;



public class MethylDbToMultisampleDomains {

	public int STEP = (int)1E6;
	public int MINCOORD = 0;
	public int MAXCOORD = (int)2.8E8;
	
	
	private static final String C_USAGE = "Use: MethylDbToDomains -outPrefix out -table methylCGsRich_normal010310_ -table methylCGsRich_tumor011010_ " + 
		" -windSize 100000 -minCpgs 6 -minMeth 0.4 -maxMeth 0.6 -minCTreads 2 -maxOppStrandAfrac 0.10 -noNonconvFilter";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-debug",usage="only does a test segment, more debug output")
    protected boolean debug = false;
    @Option(name="-debugDomain",usage="only does a test segment")
    protected boolean debugDomain = false;
    @Option(name="-table",usage="Prefix for DB table (default methylCGsRich_normal010310_)")
    protected List<String> tables = new ArrayList<String>();
    @Option(name="-withinFeat",usage="A featType from the features table")
    protected String withinFeat = null;
    @Option(name="-outPrefix",usage="Output files will have this name")
    protected String outPrefix = "wiggleTester";
    @Option(name="-minOutputWindSize",usage="only output windows this big or bigger (0)")
    protected int minOutputWindSize = 0;
    @Option(name="-windSize",usage="starting window size (500)")
    protected int windSize = 500;
    @Option(name="-variableWindowMaxWind",usage="If set , we expand window up to this size to get at least minCpgs (required but slower)")
    protected int variableWindowMaxWind = -1;
    @Option(name="-minCpgs",usage="minimum number of Cpgs in window")
    protected int minCpgs = 8;
    @Option(name="-minOutputReads",usage="don't output CpGs with less than this many reads")
    protected int minOutputReads = 0;
    @Option(name="-minMeth",usage="minimum average methylation for domain (0.4)")
    protected double minMeth = -1.0;
    @Option(name="-maxMeth",usage="maximum average methylation for domain (0.6)")
    protected double maxMeth = -1.0;
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
		new MethylDbToMultisampleDomains().doMain(args);
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
				parser.printUsage(System.err);
				System.exit(1);
			}

			if( arguments.size() != 0 ) {
				System.err.println(C_USAGE);
				parser.printUsage(System.err);
				System.exit(1);
			}

			if( (this.minMeth<0) || (this.maxMeth<0) || (this.minMeth>this.maxMeth))
			{
				System.err.println("Must provide a -minMeth and -maxMeth, and maxMeth must be larger");
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
		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.SEVERE);
		
		int nTables = this.tables.size();
		
		

		// We do fixed step or fixed wind
		// depending on input params
		CpgWalkerParams walkerParams = new CpgWalkerParams();
		walkerParams.debug = this.debug;
		walkerParams.minScanningWindCpgs = this.minCpgs;
		walkerParams.minOutputWindSize = this.minOutputWindSize;
		walkerParams.minReadsForOutput = this.minOutputReads;
		if (this.variableWindowMaxWind >= 0)
		{
			if (this.windSize>this.variableWindowMaxWind)
			{
				System.err.println("-windSize must be smaller than -variableWindowMaxWind");
				System.err.println(C_USAGE);
				parser.printUsage(System.err);
				System.err.println();
				return;
			}
				
			walkerParams.minScanningWindSize = this.windSize;
			walkerParams.maxScanningWindSize = this.variableWindowMaxWind;
			walkerParams.useVariableWindow = true;
		}
		else
		{
			walkerParams.maxScanningWindSize = this.windSize;
		}
		
		// Setup output file
		List<PrintWriter> pws = new ArrayList<PrintWriter>();
		String fixedSec = (walkerParams.useVariableWindow) ? String.format(".varMaxWind%d", variableWindowMaxWind) : "";
		String outFn = String.format("%s%s.wind%d.minOutput%d.minCpg%d.minOutReads%d.meth%.2f-%.2f.gtf", 
				this.outPrefix, fixedSec, this.windSize, this.minOutputWindSize, this.minCpgs, this.minOutputReads, this.minMeth, this.maxMeth);
		PrintWriter pw = new PrintWriter(new FileOutputStream(outFn));
		pws.add(pw);

		// Header
		pw.printf("track name=\"%s\" description=\"%s\" useScore=0 itemRgb=On visibility=4\n",
				outFn, outFn);
		
		
			
		// And domain finder
		CpgWalkerMultisampleMethRange_DEPR domainFinder = 
			new CpgWalkerMultisampleMethRange_DEPR(walkerParams,nTables, pw, this.minMeth, this.maxMeth);

		
		MethylDbQuerier params = new MethylDbQuerier();
		params.setMinCTreads(this.minCTreads);
		params.setUseNonconversionFilter(!this.noNonconvFilter);
		params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
		params.setMaxNextNonGfrac(this.maxNextNonGfrac);
		if (this.withinFeat!=null) params.addFeatFilter(this.withinFeat, 0);



		List<String> chrs = (this.debug) ? Arrays.asList("chr11","chr12") : MethylDbUtils.CHROMS;
		for (String chr : chrs)
		{
			
			
			// Set chromosomes on the domain finders
			domainFinder.setCurChr(chr);
		

			// Iterator uses DB connection and can use a ton of memory because
			// it loads all rows at once.  This stuff should really be added to iterator
			// class, but until it is , just iterate here over the chromosome
			int onCpg = 0;
			
			if (this.debugDomain)
			{
				MINCOORD = 31200000;
				MAXCOORD = 32600000;
				MINCOORD = 200000;
				MAXCOORD = 2000000;
				STEP = Math.min(STEP, MAXCOORD - MINCOORD + 1);
			}
			
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
					domainFinder.streamCpg(cpgs);

					//
					if ((onCpg % 1E5)==0) System.err.printf("On Cpg #%d, domain %s\n", onCpg, domainFinder.windStr());
					onCpg++;
				}


			}
			
			// Finish chromosomes on the domain finders
			domainFinder.finishChr();

			
		}

			
		for (PrintWriter pwi : pws)
		{
			pwi.close();
		}
		
	} // Main

	public static double safeVal(double v)
	{
		return ((Double.isNaN(v) || Double.isInfinite(v)) ? 0.0 : v);
	}
	
}


