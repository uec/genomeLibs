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
import edu.usc.epigenome.genomeLibs.MethylDb.MethylReadCollection;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylReadCollectionTiler;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairs;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsAutocorrByread;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsBinnedAutocorr;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsPearsonAutocorr;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;



public class MethylDbToTiledReads {

	public final static int STEP = (int)1E6;
	
	
	
	private static final String C_USAGE = "Use: MethylDbToTiledReads -outPrefix tiledReadTester -table methylReadCG_normalM030510_ " + 
		" -minCTreads 2 -maxOppStrandAfrac 0.10 -maxNextNonGfrac 0.10 -noNonconvFilter chr11 10000 20000";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-table",usage="Prefix for DB table (default null)")
    protected String table = null;
    @Option(name="-withinFeat",usage="A featType from the features table")
    protected String withinFeat = null;
    @Option(name="-featFlank",usage="Flank size to use with the feature")
    protected int featFlank = 0;
    @Option(name="-outPrefix",usage="Output files will have this name")
    protected String outPrefix = "tiledReadTester";
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
	
	
	protected String targetChr = null;
	protected int targetStart = -1;
	protected int targetEnd = -1;
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToTiledReads().doMain(args);
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

			if( arguments.size() != 3 ) {
				System.err.println(C_USAGE);
				System.exit(1);
			}

			targetChr = arguments.get(0);
			targetStart = Integer.parseInt(arguments.get(1));
			targetEnd = Integer.parseInt(arguments.get(2));

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
		


		// Setup streaming params
		MethylDbQuerier params = new MethylDbQuerier();
		params.setMethylTablePrefix(this.table);
		params.setReadTableSchema();
		params.setMinCTreads(this.minCTreads);
		params.setUseNonconversionFilter(!this.noNonconvFilter);
		params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
		params.setMaxNextNonGfrac(this.maxNextNonGfrac);


		// For the background, it makes a big difference whether we add this filter
		// before or after doing the 1st pass to calculate mean/sd.
		if (withinFeat!=null) params.addFeatFilter(withinFeat, this.featFlank);




		String featSec = String.format(".%s-f%d", (withinFeat==null)?"all":withinFeat, this.featFlank);
		String outFn = String.format("%s%s.%s.minCT%d.maxNextNonG%d.maxOppStrandA%d.svg",
				this.outPrefix, featSec,this.table, this.minCTreads, Math.round(100.0*this.maxNextNonGfrac), Math.round(100.0*this.maxOppStrandAfrac));


		PrintWriter pw = new PrintWriter(new FileOutputStream(outFn));

		MethylReadCollectionTiler tiler = getTiler(params, this.targetChr, this.targetStart, this.targetEnd);
		tiler.writeTiling(pw);
		
		pw.close();
		
	} // Main

	
	/**
	 * @param params
	 * @param walkerParams
	 * @param autocorrs If passed in , we use these ones.  If null, we create Pearson autocorrs from scratch,
	 *        and run them using mean 0 SD 0 (useful for doing a first pass to calculate mean/sd for second pass)
	 * @return
	 * @throws Exception
	 */
	private static MethylReadCollectionTiler getTiler(MethylDbQuerier params, String chr, int startCoord, int endCoord)
	throws Exception
	{

		System.err.printf("About to load window:%s:%d-%d\n",chr, startCoord,endCoord);
		
		MethylReadCollectionTiler reads = new MethylReadCollectionTiler(params);

		// Iterator uses DB connection and can use a ton of memory because
		// it loads all rows at once.  This stuff should really be added to iterator
		// class, but until it is , just iterate here over the chromosome
		int onCpg = 0;

		int thisStep = Math.min(STEP, endCoord-startCoord+1);


		for (int c = startCoord; c < endCoord; c += thisStep)
		{
			System.err.printf("LOADING NEW WIND: %d-%d\n",c,c+thisStep-1);

			params.clearRangeFilters();
			params.addRangeFilter(chr, c, c+thisStep-1);
			CpgIterator cpgit = new CpgIterator(params);



			//int numCpgs = cpgit.getCurNumRows();
			while (cpgit.hasNext())
			{
				Cpg cpg = cpgit.next();
				reads.addCpg(cpg);
				if ((onCpg % 1E3)==0) System.err.printf("On Cpg #%d, pos=%d, numReads=%d\n", onCpg, cpg.chromPos, reads.numReads());
				onCpg++;
			}


		}
		
		System.err.printf("Finished Making read collection, numReads=%d\n", reads.numReads());
		return reads;
	}
	

	public static double safeVal(double v)
	{
		return ((Double.isNaN(v) || Double.isInfinite(v)) ? 0.0 : v);
	}
	
}


