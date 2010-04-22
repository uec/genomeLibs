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
import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerDomainFinder;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerDomainFinderMethRange;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;



public class MethylDbToNucleosomeComparisons {

	public int STEP = (int)1E6;

	
	public int MINCOORD = 0;
	public int MAXCOORD = (int)2.8E8;
	
	public int PERIODICITY = 185;
	
	public final double LOG_OF_TWO = Math.log(2.0);
	
	public static final String methPrefix = "methylCGsRich_imr90_";
	public static final String controlMethPrefix = "methylCGsRich_normal010310_";
	
	private static final String C_USAGE = "Use: MethylDbToNucleosomeComparisons -outPrefix out  " + 
		" -withinFeat CTCF -footprintSize 145 -assocSize 20 -maxOppStrandAfrac 0.10 " + 
		" -maxNextNonGfrac 0.10 -noNonconvFilter -mnasePrefix nucTable1 -mnasePrefix nucTable2";
	
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
    protected List<String> mnasePrefixes = new ArrayList<String>(5); //"mnaseIMR90_021110_";
    @Option(name="-minReadsPerBp",usage="only output CpGs with this number of reads per bp or greater (default 0.0)")
    protected double minReadsPerBp = 0.0;
    @Option(name="-normalizationWindow",usage="window for normalWind output column (default 500)")
    protected int normalizationWindow = 500;
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
	
	public static void main(String[] args)
	throws Exception	
	{
		new MethylDbToNucleosomeComparisons().doMain(args);
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

			if( (outPrefix == null) || (this.mnasePrefixes.size()==0) ) {
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
		String mnasePrefixStr = ListUtils.excelLine(mnasePrefixes);
		String normString = (this.normalizationWindow==0) ? "" : String.format(".normalized%dbpWind", this.normalizationWindow); 
				
				
		String name = String.format("nucleosomeReads.%s.nuc%d.assoc%d%s.%s.%s",
				this.outPrefix, this.footprintSize, this.assocSize, normString, mnasePrefixStr, methPrefix);
		String outFn = String.format("%s.nuc%d.assoc%d%s.%s.%s.csv", 
				this.outPrefix, this.footprintSize, this.assocSize, normString, mnasePrefixStr, methPrefix);
		String outFnWig = outFn.replace(".csv", ".wig");
		PrintWriter pw = new PrintWriter(new FileOutputStream(outFn));
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

		int nSeries = mnasePrefixes.size();
		boolean useChipseq[] = new boolean[nSeries];
		for (int i = 0; i < nSeries; i++)
		{
			// Wow talk about a special case.
			useChipseq[i] = !(this.mnasePrefixes.get(i).contains("mnase"));
		}
		
		
		MnaseOutput readCounts[] = new MnaseOutput[nSeries];
		
		for (String chr : MethylDbUtils.CHROMS) //Arrays.asList("chr10")) // 
		{
			String s = String.format("variableStep\tchrom=%s\n", chr);
			pwWig.append(s);
			int chrInt = (new ChromFeatures()).chrom_from_public_str(chr);
		
			// Iterator uses DB connection and can use a ton of memory because
			// it loads all rows at once.  This stuff should really be added to iterator
			// class, but until it is , just iterate here over the chromosome
			int onCpg = 0;

//			MINCOORD = (int)3.5E7;
//			MAXCOORD = (int)9.7E7;
//			STEP = Math.min(STEP, MAXCOORD-MINCOORD+1);
			
			for (int c = MINCOORD; c < MAXCOORD; c += STEP)
			{
				System.err.printf("LOADING NEW WIND: %d-%d\n",c,c+STEP-1);

				params.clearRangeFilters();
				params.addRangeFilter(chr, c, c+STEP-1); 

				try
				{
				
				CpgIteratorMultisample cpgit = new CpgIteratorMultisample(params, methTables);
				//int numCpgs = cpgit.getCurNumRows();
				while (cpgit.hasNext())
				{
					Cpg[] cpgs = cpgit.next();

					Cpg cpg = cpgs[0];
					double meth = cpg.fracMeth(!this.noNonconvFilter);
	
					boolean enoughReads = false;
					for (int i = 0; i < nSeries; i++)
					{
						
						MnaseOutput mnaseReads = this.countMnaseReads(chr, cpg, mnasePrefixes.get(i), useChipseq[i]);
						// if (mnaseReads >= this.minReadsPerBp) enoughReads = true;
						enoughReads = true;
						readCounts[i] = mnaseReads;
					}
					
					//
					if ((onCpg % 1E5)==0) System.err.printf("On Cpg #%d, meth=%f, mnase=%.3f\n", onCpg, meth, readCounts[0].rawCounts);
					onCpg++;
					
					if (enoughReads)
					{
						pwWig.printf("%d\t%d\n", cpg.chromPos, Math.round(readCounts[0].val));
						pw.printf("%d,%d,%.3f,%.1f", chrInt, cpg.chromPos, meth*100,cpg.getCpgWeight());
						for (int i = 0; i < nSeries; i++)
						{
							pw.printf(",%.3f", readCounts[i].val);
						}
						pw.printf(",%d", cpgs[1].totalReads);
						for (int i = 0; i < nSeries; i++)
						{
							pw.printf(",%.3f", readCounts[i].rawCounts);
						}
						for (int i = 0; i < nSeries; i++)
						{
							pw.printf(",%.3f", readCounts[i].normWindCount);
						}
						pw.println();
					}
				}

				}
				catch (Exception e)
				{
					System.err.printf("%s\nCouldn't do region %s:%d-%d\n",e.toString(),chr,c,c+STEP-1);
					e.printStackTrace();
				}
			}
		}
		

		pw.close();
		pwWig.close();

	} // Main


	protected MnaseOutput countMnaseReads(String chr, Cpg target, String mnasePrefix, boolean chipSeqMetric)
	throws Exception
	{
		double PSEUDOCOUNT = 0.005; // per bp
		int FRAGLENHIGH = 600;
		
		
		int center = target.chromPos;
		int quarterNuc = (int)((double)this.footprintSize/4.0);
		int halfNuc = (int)((double)this.footprintSize/2.0);
		int assocHalf = (int)((double)this.assocSize/2.0);
		int fraghalf = (int)((double)FRAGLENHIGH/2.0);

		int cycle = halfNuc;
		
		int s, e;
		
		MnaseOutput out = new MnaseOutput();

		
		int windHalf = (int)((double)this.normalizationWindow/2.0);
		s = center - windHalf;
		e = center + windHalf;
		out.normWindCount = (double)this.countUnstrandedReads(mnasePrefix, chr, s, e);
		//			System.err.printf("NormWind Cpg %d\t Getting (all) reads %s at %d-%d\t%d\n", center, mnasePrefix, s, e,(int)normWindCount);

		
		if (chipSeqMetric)
		{
			// CHIP-SEQ .  THIS SHOULD TAKE  A FRAGMENT LEN, BUT 600 is based on K4 IMR90 offsets.
			
			// First fw
			s = center - fraghalf;
			e = center;
			int preCountPos = countStrandedReads(mnasePrefix, chr, s, e, StrandedFeature.POSITIVE);
			//if (preCountPos>0) System.err.printf("Cpg %d\t Getting (+) reads at %d-%d\t%d\n", center, s, e,preCountPos);

			// Then rev
			s = center;
			e = center + fraghalf;
			int postCountNeg = countStrandedReads(mnasePrefix, chr, s, e, StrandedFeature.NEGATIVE);
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
			int preCountPos = countStrandedReads(mnasePrefix, chr, s, e, StrandedFeature.POSITIVE);
//			System.err.printf("Cpg %d\t Getting (+) reads at %d-%d\t%d\n", center, s, e,preCountPos);

			// Then rev
			s = center + cycle - assocHalf;
			e = center + cycle + assocHalf;
			int postCountNeg = countStrandedReads(mnasePrefix, chr, s, e, StrandedFeature.NEGATIVE);
//			System.err.printf("Cpg %d\t Getting (-) reads at %d-%d\t%d\n", center, s, e,postCountNeg);

			int phasedCount = preCountPos + postCountNeg;
			int phasedLen = assocHalf;


			out.rawCounts = phasedCount;
			// Then unphased
			int unphasedCount = 0;
			s = center - assocHalf;
			e = center + assocHalf;
			unphasedCount += countUnstrandedReads(mnasePrefix, chr, s, e);
			s = center - this.PERIODICITY - assocHalf;
			e = center - this.PERIODICITY + assocHalf;
			unphasedCount += countUnstrandedReads(mnasePrefix, chr, s, e);
			s = center + this.PERIODICITY - assocHalf;
			e = center + this.PERIODICITY + assocHalf;
			unphasedCount += countUnstrandedReads(mnasePrefix, chr, s, e);
			int unphasedLen = 6 * assocHalf;
			//			System.err.printf("Cpg %d\t Getting unphased reads at %d-%d\t%d\n", center, center - this.PERIODICITY - assocHalf,
			//					center + this.PERIODICITY + assocHalf,unphasedCount);



			out.val = ((double)(phasedCount)+(PSEUDOCOUNT*(double)phasedLen)) * 
			(1.0/((double)unphasedCount+(PSEUDOCOUNT*(double)unphasedLen))) *
			(((double)unphasedLen)/((double)phasedLen));
			out.val = Math.log(out.val)/LOG_OF_TWO;
			//System.err.printf("\tRatio=%.1f\n", out);

//				if (this.normalizationWindow==0.0)
//				{
//
//				}
//				else
//
//				{
//					out.val = ((double)(phasedCount)+(PSEUDOCOUNT*(double)phasedLen)) * 
//					(1.0/((double)normWindCount+(PSEUDOCOUNT*(double)this.normalizationWindow))) *
//					(((double)this.normalizationWindow)/((double)phasedLen));
//					out.val = Math.log(out.val)/LOG_OF_TWO;
//					//System.err.printf("\tRatio=%.1f\n", out);
//				}

		}
		
		
		return out;
	}
	
	protected int countUnstrandedReads(String mnasePrefix, String chr, int s, int e)
	throws Exception
	{
		return countStrandedReads(mnasePrefix, chr, s, e, StrandedFeature.UNKNOWN);
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
	protected int countStrandedReads(String mnasePrefix, String chr, int s, int e, StrandedFeature.Strand targetStrand)
	throws Exception
	{
		MethylDbQuerier params = new MethylDbQuerier();
		params.setMinCTreads(1);
		params.setUseNonconversionFilter(true);
		params.setMaxOppstrandAfrac(1.0);
		params.setMaxNextNonGfrac(1.0);
		params.setMethylTablePrefix(mnasePrefix);

		params.clearRangeFilters();
		params.addRangeFilter(chr, s, e);
		CpgIterator cpgit = new CpgIterator(params);

		int count = 0;
		while (cpgit.hasNext())
		{
			// Note if you use a methylCg table, you will get CG density
			Cpg cpg = cpgit.next();
			if ((targetStrand == StrandedFeature.UNKNOWN) ||
					(cpg.getStrand() == targetStrand)) count++;
		}
		return count;
	}
	
	class MnaseOutput
	{
		public double val = 0.0;
		public double rawCounts = 0.0;
		public double normWindCount = 0.0;
	}
	
}


