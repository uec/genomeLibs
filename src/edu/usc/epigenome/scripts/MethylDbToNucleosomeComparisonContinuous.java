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

public class MethylDbToNucleosomeComparisonContinuous {

//	public int STEP = (int)1E6;
//
//	
//	public int MINCOORD = 0;
//	public int MAXCOORD = (int)2.8E8;
	
	public int PERIODICITY = 185;
	public final double LOG_OF_TWO = Math.log(2.0);
	
	public static final String controlMethPrefix = "methylCGsRich_normal010310_";
	
	private static final String C_USAGE = "Use: MethylDbToNucleosomeComparisons -outPrefix out  " + 
		" -withinFeat CTCF -footprintSize 145 -assocSize 20 -maxOppStrandAfrac 0.10 " + 
		" -maxNextNonGfrac 0.10 -noNonconvFilter -mnasePrefix nucTable1 -step 50 ";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-withinFeat",usage="A featType from the features table")
    protected String withinFeat = null;
    @Option(name="-featFlank",usage="If withinFeat is selected, we use this much flank around each feature (0)")
    protected int featFlank = 0;
    @Option(name="-outPrefix",usage="Output files will have this name")
    protected String outPrefix = null;
    @Option(name="-footprintSize",usage="nucleosome size (185)")
    protected int footprintSize = 185;
    @Option(name="-mnasePrefix",usage="mnase database table (default mnaseIMR90_021110_)")
//  protected List<String> mnasePrefixes = new ArrayList<String>(5); //"mnaseIMR90_021110_";
    protected String mnasePrefix = "mnaseIMR90_021110_";
    @Option(name="-methPrefix",usage="methylation database table (default methylCGsRich_imr90_)")
    protected String methPrefix = "methylCGsRich_imr90_";
    @Option(name="-motif",usage="If set, we use counts for this motif rather than methylation as the feature of interest (For instance 'CG' for CG density")
    protected String motif = null;
    @Option(name="-autoMnase",usage="If set, we use counts for the MNase reads themselves rather than methylation")
    protected boolean autoMnase = false;
   @Option(name="-autoMnaseFw",usage="If set, we use counts for the MNase reads themselves rather than methylation")
    protected boolean autoMnaseFw = false;
    @Option(name="-autoMnaseRev",usage="If set, we use counts for the MNase (opposite strand) reads themselves rather than methylation")
    protected boolean autoMnaseRev = false;
    @Option(name="-refGenome",usage="reference genome (default hg18)")
    protected String refGenome = "hg18";
    @Option(name="-minReadsPerBp",usage="only output CpGs with this number of reads per bp or greater (default 0.0)")
    protected double minReadsPerBp = 0.0;
    @Option(name="-normalizationWindow",usage="window for normalWind output column (default 0)")
    protected int normalizationWindow = 0;
    @Option(name="-assocSize",usage="Number of bases around 1/2 nuc to count (should be less than 1/2 nuc, default 40)")
    protected int assocSize = 20;
    @Option(name="-step",usage="Step size across chromosome (50)")
    protected int step = 50;
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

	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToNucleosomeComparisonContinuous().doMain(args);
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

			if( (arguments.size() != 0) ) {
				System.err.println(C_USAGE);
				parser.printUsage(System.err);
				System.exit(1);
			}

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
		String mnasePrefixStr = mnasePrefix;
		String normString = (this.normalizationWindow==0) ? "" : String.format(".normalized%dbpWind", this.normalizationWindow); 
		String motifString = (this.motif == null) ? "" : String.format(".motif%s", this.motif.toUpperCase());
		String withinStr = (this.withinFeat == null) ? "" : String.format(".withinFeat-%s", this.withinFeat);
				
				
		String name = String.format("nucleosomeReads.%s%s%s.minReads%.2f.nuc%d.assoc%d%s.%s.%s",
				this.outPrefix, motifString, withinStr, this.minReadsPerBp, this.footprintSize, this.assocSize, normString, mnasePrefixStr, methPrefix);
		String outFn = String.format("%s%s%s.minReads%.2f.nuc%d.assoc%d%s.%s.%s.csv", 
				this.outPrefix, motifString, withinStr, this.minReadsPerBp, this.footprintSize, this.assocSize, normString, mnasePrefixStr, methPrefix);
		String outFnWig = outFn.replace(".csv", ".wig");
//		PrintWriter pw = new PrintWriter(new FileOutputStream(outFn));
		PrintWriter pwWig = new PrintWriter(new FileOutputStream(outFnWig));
		//pw.printf("%s,%s\n","Meth", "ReadsPerBp");
		pwWig.printf("track type=wiggle_0 name=\"%s\" description=\"%s\"\n",name,name);	

		
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
		boolean useChipseq[] = new boolean[nSeries];
		for (int i = 0; i < nSeries; i++)
		{
			// Wow talk about a special case.
			useChipseq[i] = !mnasePrefix.contains("mnase") && !mnasePrefix.contains("Schones");
		}
		
		double[] methCounts = new double[METHCOUNTER_LEN];
		double[] methTotals = new double[METHCOUNTER_LEN];
		
		for (String chr :  MethylDbUtils.CHROMS) //Arrays.asList("chr11")) //MethylDbUtils.SMALL_CHROMS) //  
		{
			String s = String.format("variableStep\tchrom=%s\n", chr);
			pwWig.append(s);
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
					
					int offs = 0; // 20000000; //0;
					int offe = 0; // 81000000; //0;
					System.err.printf("offs=%d, offe=%d\n",offs,offe);
					
					// The mnase counts are the same for both cases
					counts = MethylDbUtils.chromScoresReadCounts(params, chr, this.mnasePrefix, this.refGenome, offs, offe);

					
					// The meth differ if it's a motif
					if (this.autoMnase || this.autoMnaseFw || this.autoMnaseRev)
					{
						meths = counts;
					}
					else if (this.motif == null)
					{
						meths = MethylDbUtils.chromScoresMethLevels(params, chr, this.methPrefix, this.refGenome, offs, offe);
					}
					else
					{
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
						System.err.printf("About to populate both strands for %s for motif %s\n", chr, this.motif);
						all.populate(chr, this.motif, seqArr, offs, StrandedFeature.UNKNOWN);

						meths[0] = all;
						meths[1] = all;

//						System.err.printf("Fw score at %d = %d\n",24009253,counts[0].getScore(chr, 24009253));
//						WigOptions wo = new WigOptions();
//						wo.f_step = 1;
//						//counts[0] = counts[0].smooth(50, 30);
//						counts[0].wigOutput("testFw.wig", wo);
//						counts[1].wigOutput("testRev.wig", wo);

					}

					
					int minPos = counts[0].chromMinPos(chr); 
					System.err.println("Getting min pos: " + minPos);
					int maxPos = counts[0].chromMaxPos(chr);
					System.err.println("Getting max pos: " + maxPos);
					
					for (int pos = minPos; pos<maxPos; pos += 1) // this.step)
					{

						boolean fwRead = (counts[0].getScore(chr, pos).intValue() >= 1);
						boolean revRead = (counts[1].getScore(chr, pos).intValue() >= 1);
						
						boolean enoughReads = true;
						MnaseOutput mnaseReads = new MnaseOutput();
						mnaseReads.rawCounts = (fwRead) ? 1 : 0;
						if (this.minReadsPerBp > 0.0)
						{
							if (fwRead)
							{
								int nucCenter = pos + ( (this.PERIODICITY-this.assocSize) / 2);
								mnaseReads = this.countMnaseReads(chr, nucCenter, counts[0], counts[1], useChipseq[0]);
							}
							else if (revRead)
							{
								int nucCenter = pos - ( (this.PERIODICITY-this.assocSize) / 2);
								mnaseReads = this.countMnaseReads(chr, nucCenter, counts[0], counts[1], useChipseq[0]);
							}
							enoughReads = (mnaseReads.val >= this.minReadsPerBp);
						}

						
						
						//if (mnaseReads.rawCounts>2) System.err.printf("Raw counts=%d\n", (int)mnaseReads.rawCounts);


						//enoughReads = true;

						if (enoughReads && fwRead) incrementMethCounts(chr,pos,false,methCounts,methTotals,meths,maxPos);
						if (enoughReads && revRead) incrementMethCounts(chr,pos,true,methCounts,methTotals,meths,maxPos);

						
					
						//
						if ((pos % 1E6)==0) System.err.printf("On pos #%d, meth=%d\n", pos, (int)mnaseReads.rawCounts);
					
						if (enoughReads)
						{
//							pwWig.printf("%d\t%.2f\n",pos, mnaseReads.val);
////							pwWig.printf("%d\t%.2f\t%.2f\t%s\n",pos,mnaseReads.rawCounts, mnaseReads.val, fwRead?"+":"-");
//
////							pw.printf("%d,%d", chrInt, pos);
//////							pw.printf(",%.3f", mnaseReads.val);
////							pw.printf(",%d", (int)mnaseReads.rawCounts);
//////							pw.printf(",%.3f", mnaseReads.normWindCount);
////							pw.println();
						}
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
		pwWig.close();
		
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
				if (this.autoMnaseRev || this.autoMnaseFw)
				{
					methCounts[i] += 1.0;
				}
				else if (this.autoMnase)
				{
					methCounts[i] += 2.0;
				}
				else if (this.motif != null)
				{
					methCounts[i] += 2.0;
				}
				else if (cpg)
				{
					methCounts[i] += 1.0;
				}
			
				// And the score
				if (this.autoMnaseFw || this.autoMnaseRev)
				{
					int ind = (rev) ? 1 : 0;
					if (this.autoMnaseRev) ind = 1-ind;
					double m = meths[ind].getScore(chr, coord).doubleValue();
					methTotals[i] += m;
				}
				else if (this.autoMnase)
				{
					double m = meths[0].getScore(chr, coord).doubleValue() + meths[1].getScore(chr, coord).doubleValue();
					methTotals[i] += m;					
				}
				else if (cpg)
				{
					double m = meths[1].getScore(chr, coord).doubleValue();
					//if (m>1.0) System.err.printf("Meth>1.0 (%.3f)\n",m);
					methTotals[i] += m; 
				}
			}
		}
		
	}
	

	protected MnaseOutput countMnaseReads(String chr, int pos, ChromScoresFast fwCounts, ChromScoresFast revCounts, boolean chipSeqMetric)
	throws Exception
	{
		double PSEUDOCOUNT = 0.005; // per bp
		int FRAGLENHIGH = 600;
		
		
		int center = pos;
		int quarterNuc = (int)((double)this.footprintSize/4.0);
		int halfNuc = (int)((double)this.footprintSize/2.0);
		int assocHalf = (int)((double)this.assocSize/2.0);
		int fraghalf = (int)((double)FRAGLENHIGH/2.0);

		int cycle = halfNuc;
		
		int s, e;
		
		MnaseOutput out = new MnaseOutput();

		
//		int windHalf = (int)((double)this.normalizationWindow/2.0);
//		s = center - windHalf;
//		e = center + windHalf;
//		out.normWindCount = (double)this.countUnstrandedReads(fwCounts, revCounts,chr, s, e);
//		//			System.err.printf("NormWind Cpg %d\t Getting (all) reads %s at %d-%d\t%d\n", center, mnasePrefix, s, e,(int)normWindCount);

		
		if (chipSeqMetric)
		{
			// CHIP-SEQ .  THIS SHOULD TAKE  A FRAGMENT LEN, BUT 600 is based on K4 IMR90 offsets.
			System.err.println("Why are we using ChIP-seq metric??");
			System.exit(1);
			
			// First fw
			s = center - fraghalf;
			e = center;
			int preCountPos = countReads(fwCounts, chr, s, e);
			//if (preCountPos>0) System.err.printf("Cpg %d\t Getting (+) reads at %d-%d\t%d\n", center, s, e,preCountPos);

			// Then rev
			s = center;
			e = center + fraghalf;
			int postCountNeg = countReads(revCounts, chr, s, e);
			//if (postCountNeg>0) System.err.printf("Cpg %d\t Getting (-) reads at %d-%d\t%d\n", center, s, e,postCountNeg);

			out.val = (double)(preCountPos+postCountNeg)/(double)fraghalf;
			out.rawCounts = preCountPos+postCountNeg;
			
			
//			if (this.normalizationWindow>0.0)
//			{
//				out = ((double)(preCountPos+postCountNeg)+(PSEUDOCOUNT*(double)fraghalf)) * 
//					(1.0/((double)normWindCount+(PSEUDOCOUNT*(double)this.normalizationWindow))) *
//							(((double)this.normalizationWindow)/((double)fraghalf));
//				out = Math.log(out)/LOG_OF_TWO;
//			}
		}
		else
		{
			// MNASE METRIC

			// First fw
			s = center - cycle - assocHalf;
			e = center - cycle + assocHalf;
			int preCountPos = countReads(fwCounts, chr, s, e);
			//System.err.printf("Cpg %d\t Getting (+) reads at %d-%d\t%d\n", center, s, e,preCountPos);

			// Then rev
			s = center + cycle - assocHalf;
			e = center + cycle + assocHalf;
			int postCountNeg = countReads(revCounts, chr, s, e);
			//System.err.printf("Cpg %d\t Getting (-) reads at %d-%d\t%d\n", center, s, e,postCountNeg);

			int phasedCount = preCountPos + postCountNeg;
			int phasedLen = 2*assocHalf;
			
			out.rawCounts = phasedCount;
			out.val = phasedCount;

			
			
			
			// Then unphased
			int unphasedCount = 0;
			s = center - assocHalf;
			e = center + assocHalf;
			unphasedCount += countUnstrandedReads(fwCounts, revCounts, chr, s, e);
			int unphasedLen = 2 * assocHalf;
//			s = center - this.PERIODICITY - assocHalf;
//			e = center - this.PERIODICITY + assocHalf;
//			unphasedCount += countUnstrandedReads(fwCounts, revCounts, chr, s, e);
//			s = center + this.PERIODICITY - assocHalf;
//			e = center + this.PERIODICITY + assocHalf;
//			unphasedCount += countUnstrandedReads(fwCounts, revCounts,  chr, s, e);
//			int unphasedLen = 6 * assocHalf;
//			//			System.err.printf("Cpg %d\t Getting unphased reads at %d-%d\t%d\n", center, center - this.PERIODICITY - assocHalf,
//			//					center + this.PERIODICITY + assocHalf,unphasedCount);
//
//
//

			out.val = ((double)(phasedCount)+(PSEUDOCOUNT*(double)phasedLen)) * 
			(1.0/((double)unphasedCount+(PSEUDOCOUNT*(double)unphasedLen))) *
			(((double)unphasedLen)/((double)phasedLen));
			out.val = Math.log(out.val)/LOG_OF_TWO;

			out.val = phasedCount - unphasedCount;
			out.val = Math.min(preCountPos, postCountNeg) - unphasedCount;


		}
		
		
		return out;
	}
	
	protected int countUnstrandedReads(ChromScoresFast fw, ChromScoresFast rev, String chr, int s, int e)
	throws Exception
	{
		int out = countReads(fw,chr,s,e) + countReads(rev,chr,s,e);
		return out;
	}
	
	/**
	 * @param mnasePrefix
	 * @param chr
	 * @param s
	 * @param e
	 * @param targetStrand: Set to UNKNOWN if you want either strand
	 * @return
	 * @throws Exception
	 */
	protected int countReads(ChromScoresFast scores, String chr, int s, int e)
	throws Exception
	{
		int out = (int) scores.getScoresTotal(chr,s,e);
		//System.err.printf("Getting stranded reads %s:%d-%d = %d\n",chr,s,e,out);
		return out;
	}
	
	class MnaseOutput
	{
		public double val = 0.0;
		public double rawCounts = 0.0;
		public double normWindCount = 0.0;
	}
	
}


