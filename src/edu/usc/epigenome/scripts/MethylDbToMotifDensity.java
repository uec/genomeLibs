package edu.usc.epigenome.scripts;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresFast;
import edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresMotifPositions;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairs;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsBinnedAutocorr;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerAllpairsPearsonAutocorr;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;



public class MethylDbToMotifDensity {

	public int STEP = (int)1E6;
	public int MINCOORD = 0;
	public int MAXCOORD = (int)2.8E8;
	
	
	
	private static final String C_USAGE = "Use: MethylDbToAutocorr -outPrefix out -table methylCGsRich_normal010310_ -table methylCGsRich_tumor011010_ " + 
		" -windSize 500 -minCTreads 2 -maxOppStrandAfrac 0.10 -maxNextNonGfrac 0.10 -refGenome hg18 -noNonconvFilter AAAA";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-table",usage="Prefix for DB table (default methylCGsRich_normal010310_)")
    protected List<String> tables = new ArrayList<String>();
    @Option(name="-withinFeat",usage="A featType from the features table")
    protected List<String> withinFeats = new ArrayList<String>(10);
    @Option(name="-featFlank",usage="Flank size to use with the feature")
    protected int featFlank = 0;
    @Option(name="-outPrefix",usage="Output files will have this name")
    protected String outPrefix = "wiggleTester";
    @Option(name="-windSize",usage="motif density window size (500)")
    protected int windSize = 500;
    @Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
    protected int minCTreads = 2;
    @Option(name="-refGenome",usage="reference genome (default hg18)")
    protected String refGenome = "hg18";
    @Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
    		"This sets a maximum number of observed As on the opposite strand (default 0.1)")
    protected double maxOppStrandAfrac = 0.1;
    @Option(name="-maxNextNonGfrac",usage="If the base following the C has more than this ratio of non-G bases, we don't count it. (default 0.1)")
    protected double maxNextNonGfrac = 0.1;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	protected int halfWind = -1;
	
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToMotifDensity().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception
	{
		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		parser.setUsageWidth(80);
		String motif = null;
		try
		{
			parser.parseArgument(args);

			if( tables.size() < 1)
			{
				System.err.println("Must specify at least one table");
				System.err.println(C_USAGE);
				System.exit(1);
			}

			if( arguments.size() != 1 ) {
				System.err.println(C_USAGE);
				System.exit(1);
			}

			motif = arguments.get(0);

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
		
		this.halfWind = (int)((double)this.windSize/2.0);
		int nTables = this.tables.size();
		

		
		
		// Loop through within feats.  If more than one is specified, throw in a null one for fun.  Also add
		// it if we haven't specified any
		if (withinFeats.size() != 1) withinFeats.add(null);

		for (String withinFeat : withinFeats)
		{

			// Setup output files
			List<PrintWriter> pws = new ArrayList<PrintWriter>();

			for (int i = 0; i < nTables; i++)
			{
				// Start table
				String tab = this.tables.get(i);

				String motifString = (motif == null) ? "" : String.format(".motif%s", motif.toUpperCase());
				String featSec = String.format(".%s-f%d", (withinFeat==null)?"all":withinFeat, this.featFlank);
				String outFn = String.format("MotifDensities.%s.%s%s%s.wind%d.csv", 
						this.outPrefix, tab, featSec, motifString, this.windSize);
				PrintWriter pw = new PrintWriter(new FileOutputStream(outFn));
				pws.add(pw);
			}

			// Setup streaming params
			MethylDbQuerier params = new MethylDbQuerier();
			params.setMinCTreads(this.minCTreads);
			params.setUseNonconversionFilter(!this.noNonconvFilter);
			params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
			params.setMaxNextNonGfrac(this.maxNextNonGfrac);

			
			// For the background, it makes a big difference whether we add this filter
			// before or after doing the 1st pass to calculate mean/sd.
			if (withinFeat!=null) params.addFeatFilter(withinFeat, this.featFlank);


			for (int i = 0; i < nTables; i++)
			{
				// Header
				PrintWriter pw = pws.get(i);
				String header = String.format("totalCT,meth,%s-freq",motif);
				if (header != null) pw.println(header);
			}

			// And run!
			this.outputCGs(params, pws, motif);



			// Finish up
			for (int i = 0; i < nTables; i++)
			{
				PrintWriter pw = pws.get(i);
				pw.close();
			}

		}

		
	} // Main

	
	/**
	 * @param params
	 * @param pws One for each table
	 * @param motif
	 * @return
	 * @throws Exception
	 */
	private void outputCGs(MethylDbQuerier params, List<PrintWriter> pws, String motif) 
	throws Exception
	{
		int nTabs = this.tables.size();
			


		for (String chr : Arrays.asList("chr2","chr3")) // MethylDbUtils.CHROMS) //  
		{

			
			// Get the motif positions
			ChromScoresMotifPositions motifPos = getChromMotifPositions(chr, motif);

			// Iterator uses DB connection and can use a ton of memory because
			// it loads all rows at once.  This stuff should really be added to iterator
			// class, but until it is , just iterate here over the chromosome
			int onCpg = 0;

//			MINCOORD = 7000000;
//			MAXCOORD = 9000000;
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
					int coord = cpgs[0].chromPos;

					int nMotifs = countMotifs(motifPos, chr, coord);
					double motifFreq = (double)nMotifs / (double)this.windSize;
					
					// Output cpg
					for (int i = 0; i < nTabs; i++)
					{
						//System.err.printf("Streaming cpg at pos %d to table %s\n",cpgs[i].chromPos,this.tables.get(i));
						Cpg cpg = cpgs[i];
						
						PrintWriter pw = pws.get(i);
						pw.printf("%d,%.3f,%.3f\n",cpg.totalReads,cpg.fracMeth(!this.noNonconvFilter),motifFreq);
					}

					//
					if ((onCpg % 1E5)==0) System.err.printf("On Cpg #%d\n", onCpg);
					onCpg++;
				}


			}
		}

	}
	
	private ChromScoresMotifPositions getChromMotifPositions(String chr, String motif)
	throws Exception
	{
		int offs = 0;
		int offe = 0;
		
		char[] seqArr = null;
		{
			Sequence seq = GoldAssembly.chromSeq(this.refGenome, chr);
			String seqStr = seq.seqString();
			if (offe>0)
			{
				seqStr = seqStr.substring(offs, offe);
			}
			seqArr = seqStr.toUpperCase().toCharArray();
		}
		
		System.err.println("Seq length=" + seqArr.length);
		
		
		ChromScoresMotifPositions all = new ChromScoresMotifPositions(this.refGenome);
		System.err.printf("About to populate both strands for %s for motif %s\n", chr, motif);
		all.populate(chr, motif, seqArr, offs, StrandedFeature.UNKNOWN);		
		
		return all;
	}
		
	protected int countMotifs(ChromScoresFast scores, String chr, int coordCenter)
	throws Exception
	{
		int s = Math.max(1, coordCenter-this.halfWind);
		int e = Math.min(GoldAssembly.chromLengthStatic(chr, this.refGenome), coordCenter+this.halfWind);
		return countMotifs(scores, chr, s, e);
	}
	
	protected int countMotifs(ChromScoresFast scores, String chr, int s, int e)
	throws Exception
	{
		int out = (int) scores.getScoresTotal(chr,s,e);
		//System.err.printf("Getting stranded reads %s:%d-%d = %d\n",chr,s,e,out);
		return out;
	}

}


