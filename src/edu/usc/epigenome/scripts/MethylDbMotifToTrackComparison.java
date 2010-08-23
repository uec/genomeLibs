package edu.usc.epigenome.scripts;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.WigOptions;
import edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresArray;
import edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresFast;
import edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresMotifPositions;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerDomainFinder;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerDomainFinderMethRange;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;

public class MethylDbMotifToTrackComparison {

//	public int STEP = (int)1E6;
//
//	
//	public int MINCOORD = 0;
//	public int MAXCOORD = (int)2.8E8;
	
	public int PERIODICITY = 185;
	public final double LOG_OF_TWO = Math.log(2.0);
	
	public static final String controlMethPrefix = "methylCGsRich_normal010310_";
	
	private static final String C_USAGE = "Use: MethylDbMotifToTrackComparison -outPrefix out  " + 
		" -withinFeat CTCF -footprintSize 145 -assocSize 20 -maxOppStrandAfrac 0.10 " + 
		" -maxNextNonGfrac 0.10 -noNonconvFilter -methPrefix methylCGsRich_imr90_ -onlyFwMotifs CG ";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-onlyFwMotifs",usage="If set, we only use motifs from the fw strand of the reference genome (default false)")
    protected boolean onlyFwMotifs = false;
    @Option(name="-withinFeat",usage="A featType from the features table")
    protected String withinFeat = null;
    @Option(name="-featFlank",usage="If withinFeat is selected, we use this much flank around each feature (0)")
    protected int featFlank = 0;
    @Option(name="-outPrefix",usage="Output files will have this name")
    protected String outPrefix = null;
    @Option(name="-footprintSize",usage="nucleosome size (185)")
    protected int footprintSize = 145;
    @Option(name="-methPrefix",usage="methylation database table (default methylCGsRich_imr90_)")
    protected String methPrefix = "methylCGsRich_imr90_";
    @Option(name="-refGenome",usage="reference genome (default hg18)")
    protected String refGenome = "hg18";
    @Option(name="-minReadsPerBp",usage="only output CpGs with this number of reads per bp or greater (default 0.0)")
    protected double minReadsPerBp = 0.0;
    @Option(name="-normalizationWindow",usage="window for normalWind output column (default 0)")
    protected int normalizationWindow = 0;
    @Option(name="-assocSize",usage="Number of bases around 1/2 nuc to count (should be less than 1/2 nuc, default 40)")
    protected int assocSize = 20;
    @Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
    protected int minCTreads = 1;
    @Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
    		"This sets a maximum number of observed As on the opposite strand (default 0.1)")
    protected double maxOppStrandAfrac = 0.1;
    @Option(name="-maxNextNonGfrac",usage="If the base following the C has more than this ratio of non-G bases, we don't count it. (default 0.1)")
    protected double maxNextNonGfrac = 0.1;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();
	
	public static final int METHCOUNTER_LEN = 2000; // 1600; // 2000; //1000;
	public static final int METHCOUNTER_HALFLEN = METHCOUNTER_LEN / 2;

	protected String motif = null;
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbMotifToTrackComparison().doMain(args);
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

			if( (arguments.size() != 1) ) {
				System.err.println(C_USAGE);
				parser.printUsage(System.err);
				System.exit(1);
			}

			this.motif = arguments.get(0);
			if( (outPrefix == null) ) {
				System.err.println("Must specify outPrefix and at least one mnasePrefix");
				parser.printUsage(System.err);
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
		
		

		// Setup output files and print domain finders
		ListUtils.setDelim("-");
		String mnasePrefixStr = "";
		String strandSec = (this.onlyFwMotifs) ? ".fwStrandOnly" : ".bothStrands";
		String normString = (this.normalizationWindow==0) ? "" : String.format(".normalized%dbpWind", this.normalizationWindow); 
		String motifString = (this.motif == null) ? "" : String.format(".motif%s", this.motif.toUpperCase());
		String withinStr = (this.withinFeat == null) ? "" : String.format(".withinFeat-%s", this.withinFeat);
				
				
		String name = String.format("nucleosomeReads.%s%s%s.minReads%.2f.nuc%d.assoc%d%s.%s.%s",
				this.outPrefix, motifString, withinStr, this.minReadsPerBp, this.footprintSize, this.assocSize, strandSec, mnasePrefixStr, methPrefix);
		String outFn = String.format("%s%s%s.minReads%.2f.nuc%d.assoc%d%s.%s.%s.csv", 
				this.outPrefix, motifString, withinStr, this.minReadsPerBp, this.footprintSize, this.assocSize, strandSec, mnasePrefixStr, methPrefix);
		
		MethylDbQuerier params = new MethylDbQuerier();
		params.setMinCTreads(this.minCTreads);
		params.setUseNonconversionFilter(!this.noNonconvFilter);
		params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
		params.setMaxNextNonGfrac(this.maxNextNonGfrac);
		if (this.withinFeat!=null) params.addFeatFilter(this.withinFeat, this.featFlank);

		// We use the control table simply to limit the valid Cpg.  This is a result of 
		// our incorrect loading of Lister 2009 tables, which contains many Cpgs which 
		// incorrectly yield a meth level of 0 at CpGs not covered in their sequencing data.
		// This was an artifact of the way they published their data, where they published
		// a list of methy-C positions without positions containing 0 mC reads, so we had
		// to add fake positions for all Cs in the genome, and we didn't know which ones
		// actually had coverage in their data.
		List<String> methTables = Arrays.asList(methPrefix, controlMethPrefix);
	//	List<String> mnaseTables = Arrays.asList(mnasePrefix);

		int nSeries = 1;
		
		double[] methCounts = new double[METHCOUNTER_LEN];
		double[] methTotals = new double[METHCOUNTER_LEN];
		
		int chromNum=1;
		for (String chr :  MethylDbUtils.CHROMS_MINUS_TWELVE) //  Arrays.asList("chr21","chr22")) //
		{
			System.err.printf("On chrom %d (%s)\n",chromNum++,chr);
			String s = String.format("variableStep\tchrom=%s\n", chr);
			int chrInt = (new ChromFeatures()).chrom_from_public_str(chr);
		
			// Iterator uses DB connection and can use a ton of memory because
			// it loads all rows at once.  This stuff should really be added to iterator
			// class, but until it is , just iterate here over the chromosome
			int onCpg = 0;

			
			// Get the full array
			
				try
				{
					{
					ChromScoresFast counts[] = new ChromScoresFast[2];
					ChromScoresFast meths[] = new ChromScoresFast[2];
					
					int offs = 0; //20000000; //0;
					int offe = 0; //21000000; //0;
					System.err.printf("offs=%d, offe=%d\n",offs,offe);
					
					counts = MethylDbUtils.chromScoresMotifCounts(chr, this.refGenome, this.motif, offs, offe,this.onlyFwMotifs);
					meths = MethylDbUtils.chromScoresMethLevels(params, chr, this.methPrefix, this.refGenome, offs, offe);


					
					int minPos = counts[0].chromMinPos(chr); 
					System.err.println("Getting min pos: " + minPos);
					int maxPos = counts[0].chromMaxPos(chr);
					System.err.println("Getting max pos: " + maxPos);
					
					
					for (int pos = minPos; pos<maxPos; pos += 1) // this.step)
					{

						boolean fwRead = (counts[0].getScore(chr, pos).intValue() >= 1);
						boolean revRead = (this.onlyFwMotifs) ? false : (counts[1].getScore(chr, pos).intValue() >= 1);
						
						boolean enoughReads = true;
						
						//if (mnaseReads.rawCounts>2) System.err.printf("Raw counts=%d\n", (int)mnaseReads.rawCounts);


						//enoughReads = true;

						if (enoughReads && fwRead) incrementMethCounts(chr,pos,false,methCounts,methTotals,meths,maxPos);
						if (enoughReads && revRead) incrementMethCounts(chr,pos,true,methCounts,methTotals,meths,maxPos);

						
					
						//
						if ((pos % 1E6)==0) System.err.printf("On pos #%d\n", pos);
					

					}

					// Try to get rid of objects
					counts[0] = null;
					counts[1] = null;
					meths[0] = null;
					meths[1] = null;
					System.gc();
					}
					
				}
				catch (Exception e)
				{
					System.err.printf("%s\nCouldn't do region %s\n",e.toString(),chr);
					e.printStackTrace();
				}
		}
	
		

//		pw.close();
		
		outFn = outFn.replace(".csv", ".methAlign.csv");
		PrintWriter pw = new PrintWriter(new FileOutputStream(outFn));
		ListUtils.setDelim(",");
		pw.println(ListUtils.excelLine(methCounts));
		pw.println(ListUtils.excelLine(methTotals));
		pw.close();

	} // Main

	/**
	 * @param chr
	 * @param pos
	 * @param rev
	 * @param methCounts
	 * @param methTotals
	 * @param meths meths[0] is the presence/absence of a CpG, meths[1] is the CpG level.
	 * @param maxCoord
	 */
	protected void incrementMethCounts(String chr, int pos, boolean rev, double[] methCounts, double[] methTotals, ChromScoresFast meths[], int maxCoord)
	{
		int len = methCounts.length;
		
		for (int i = 0; i < len; i++)
		{
			int coord = (rev) ? ((pos+METHCOUNTER_HALFLEN) - i) : ((pos-METHCOUNTER_HALFLEN) + i);
			if ((coord>=0) && (coord<maxCoord))
			{
				boolean cpg = meths[0].getScore(chr, coord).doubleValue() > 0.0;
				
				// Increment the totals count
				if (cpg) methCounts[i] += 1.0;

			
				// And the score
				if (cpg)
				{
					double m = meths[1].getScore(chr, coord).doubleValue();
					//if (m>1.0) System.err.printf("Meth>1.0 (%.3f)\n",m);
					methTotals[i] += m; 
				}
			}
		}
		
	}
	


	
}


