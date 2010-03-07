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
import edu.usc.epigenome.genomeLibs.MethylDb.CpgMethDiff;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerDomainFinder;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerDomainFinderMethRange;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;



public class MethylDbToDiffDomains {

	public int STEP = (int)1E6;
	public int MAXCOORD = (int)2.8E8;
	
	
	private static final String C_USAGE = "Use: MethylDbToDiffDomains -outPrefix out -table methylCGsRich_normal010310_ " + 
		" -windSize 100000 -minCpgs 6 -minMeth 0.4 -maxMeth 0.6 -minCTreads 2 -maxOppStrandAfrac 0.10 -noNonconvFilter table1 table2";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-withinFeat",usage="A featType from the features table")
    protected String withinFeat = null;
    @Option(name="-outPrefix",usage="Output files will have this name")
    protected String outPrefix = "wiggleTester";
    @Option(name="-windSize",usage="starting window size (500)")
    protected int windSize = 500;
    @Option(name="-minCpgs",usage="starting window size (500)")
    protected int minCpgs = 8;
    @Option(name="-minMeth",usage="minimum average methylation for domain (0.4)")
    protected double minMeth = 0.4;
    @Option(name="-maxMeth",usage="maximum average methylation for domain (0.6)")
    protected double maxMeth = 0.6;
    @Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
    protected int minCTreads = 2;
    @Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
    		"This sets a maximum number of observed As on the opposite strand (default 0.1)")
    protected double maxOppStrandAfrac = 0.1;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToDiffDomains().doMain(args);
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


			if( arguments.size() != 2 ) {
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
		List<String> tables = arguments;
//		int nTables = tables.size();
		
		

		// Setup output files and print domain finders
		List<PrintWriter> pws = new ArrayList<PrintWriter>();
		CpgWalkerParams walkerParams = new CpgWalkerParams();
		walkerParams.maxScanningWindSize = this.windSize;
		walkerParams.minScanningWindCpgs = this.minCpgs;

		String outFn = String.format("%s.DiffA-%s.DiffB-%s.wind%d.minCpg%d.meth%.2f-%.2f.bed", 
				this.outPrefix, tables.get(0), tables.get(1), this.windSize, this.minCpgs, this.minMeth, this.maxMeth);
		PrintWriter pw = new PrintWriter(new FileOutputStream(outFn));
		pws.add(pw);

		// Header
		pw.printf("track name=\"%s\" description=\"%s\" useScore=0 itemRgb=On visibility=4\n",
				outFn, outFn);
		CpgWalkerDomainFinder domainFinder = new CpgWalkerDomainFinderMethRange(walkerParams, null, pw,
				this.minMeth, this.maxMeth);
		
		// And the querier params
		MethylDbQuerier params = new MethylDbQuerier();
		params.setMinCTreads(this.minCTreads);
		params.setUseNonconversionFilter(!this.noNonconvFilter);
		params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
		if (this.withinFeat!=null) params.addFeatFilter(this.withinFeat, 0);

		for (String chr : MethylDbUtils.CHROMS) //Arrays.asList("chr1")) //  
		{
			// Set chromosomes on the domain finders
			domainFinder.setCurChr(chr);

			// Iterator uses DB connection and can use a ton of memory because
			// it loads all rows at once.  This stuff should really be added to iterator
			// class, but until it is , just iterate here over the chromosome
			int onCpg = 0;
			for (int c = 0; c < MAXCOORD; c += STEP)
//			MAXCOORD = 2981976;
//			for (int c = 2964914 ; c < MAXCOORD; c+=STEP)
			{
				System.err.printf("LOADING NEW WIND: %d-%d\n",c,c+STEP-1);

				params.clearRangeFilters();
				int last = Math.min(c+STEP-1, MAXCOORD);
				params.addRangeFilter(chr, c, last);
				CpgIteratorMultisample cpgit = new CpgIteratorMultisample(params, tables);
				//int numCpgs = cpgit.getCurNumRows();
				while (cpgit.hasNext())
				{
					Cpg[] cpgs = cpgit.next();
					
					// Make a consensus CpG where the methyl level is the differential.
					CpgMethDiff diffCpg = new CpgMethDiff(cpgs[0],cpgs[1]);

					domainFinder.streamCpg(diffCpg);

					//
					if ((onCpg % 1E5)==0) System.err.printf("On Cpg #%d, domain %s\n", onCpg, domainFinder.windStr());
					onCpg++;
				}


			}
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


