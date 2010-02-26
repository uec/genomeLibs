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
	public int MAXCOORD = (int)2.8E8;
	
	public static final String methPrefix = "methylCGsRich_imr90_";
	public static final String controlMethPrefix = "methylCGsRich_normal010310_";
	
	private static final String C_USAGE = "Use: MethylDbToNucleosomeComparisons -outPrefix out  " + 
		" -withinFeat CTCF -footprintSize 145 -assocSize 20 -maxOppStrandAfrac 0.10 " + 
		" -maxNextNonGfrac 0.10 -noNonconvFilter -mnasePrefix nucTable1 -mnasePrefix nucTable2";
	
    @Option(name="-noNonconvFilter",usage="override the nonconversion filter (default false)")
    protected boolean noNonconvFilter = false;
    @Option(name="-withinFeat",usage="A featType from the features table")
    protected String withinFeat = null;
    @Option(name="-outPrefix",usage="Output files will have this name")
    protected String outPrefix = null;
    @Option(name="-footprintSize",usage="nucleosome size (185)")
    protected int footprintSize = 185;
    @Option(name="-mnasePrefix",usage="mnase database table (default mnaseIMR90_021110_)")
    protected List<String> mnasePrefixes = new ArrayList<String>(5); //"mnaseIMR90_021110_";
    @Option(name="-minReadsPerBp",usage="only output CpGs with this number of reads per bp or greater (default 0.0)")
    protected double minReadsPerBp = 0.0;
    @Option(name="-assocSize",usage="Number of bases around 1/2 nuc to count (should be less than 1/2 nuc, default 40)")
    protected int assocSize = 40;
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
		String outFn = String.format("%s.nuc%d.assoc%d.%s.%s.csv", 
				this.outPrefix, this.footprintSize, this.assocSize, mnasePrefixStr, methPrefix);
		PrintWriter pw = new PrintWriter(new FileOutputStream(outFn));
		pw.printf("%s,%s\n","Meth", "ReadsPerBp");

		
		MethylDbQuerier params = new MethylDbQuerier();
		params.setMinCTreads(this.minCTreads);
		params.setUseNonconversionFilter(!this.noNonconvFilter);
		params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
		params.setMaxNextNonGfrac(this.maxNextNonGfrac);
		if (this.withinFeat!=null) params.addFeatFilter(this.withinFeat, 0);

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
		double readCounts[] = new double[nSeries];
		
		for (String chr : MethylDbUtils.CHROMS) //Arrays.asList("chr2","chr21","chr22")) // 
		{
			
			// Iterator uses DB connection and can use a ton of memory because
			// it loads all rows at once.  This stuff should really be added to iterator
			// class, but until it is , just iterate here over the chromosome
			int onCpg = 0;
			for (int c = 0; c < MAXCOORD; c += STEP)
//			for (int c = 0; c < 1000000; c += STEP)
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
						double mnaseReads = this.countMnaseReads(chr, cpg, mnasePrefixes.get(i));
						if (mnaseReads >= this.minReadsPerBp) enoughReads = true;
						readCounts[i] = mnaseReads;
					}
					
					//
					if ((onCpg % 1E5)==0) System.err.printf("On Cpg #%d, meth=%f, mnase=%.3f\n", onCpg, meth, readCounts[0]);
					onCpg++;
					
					if (enoughReads)
					{
						//pw.printf("%s,%d,%.3f,%.3f\n", chr, cpg.chromPos, meth*100, mnaseReads);
						pw.printf("%.3f", meth*100);
						for (int i = 0; i < nSeries; i++)
						{
							pw.printf(",%.3f", readCounts[i]);
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

	} // Main


	protected double countMnaseReads(String chr, Cpg target, String mnasePrefix)
	throws Exception
	{
		MethylDbQuerier params = new MethylDbQuerier();
		params.setMinCTreads(1);
		params.setUseNonconversionFilter(true);
		params.setMaxOppstrandAfrac(1.0);
		params.setMaxNextNonGfrac(1.0);
		params.setMethylTablePrefix(mnasePrefix);

		int center = target.chromPos;
		int quarterNuc = (int)((double)this.footprintSize/4.0);
		int halfNuc = (int)((double)this.footprintSize/2.0);
		int assocHalf = (int)((double)this.assocSize/2.0);

		int cycle = halfNuc;
		
		int s, e;
		CpgIterator cpgit;
		
		// First fw
		s = center - cycle - assocHalf;
		e = center - cycle + assocHalf;
		params.clearRangeFilters();
		params.addRangeFilter(chr, s, e);
		cpgit = new CpgIterator(params);
		int fwCount = countStrandedReads(cpgit, StrandedFeature.POSITIVE);
		//System.err.printf("Cpg %d\t Getting (+) reads at %d-%d\t%d\n", center, s, e,fwCount);

		// Then rev
		s = center + cycle - assocHalf;
		e = center + cycle + assocHalf;
		params.clearRangeFilters();
		params.addRangeFilter(chr, s, e);
		cpgit = new CpgIterator(params);
		int revCount = countStrandedReads(cpgit, StrandedFeature.NEGATIVE);
		//System.err.printf("Cpg %d\t Getting (-) reads at %d-%d\t%d\n", center, s, e,revCount);
		
		double frac = ((double)fwCount+(double)revCount)/(2.0*(double)assocHalf);
		
		return frac;
	}
	
	protected int countStrandedReads(CpgIterator cpgit, StrandedFeature.Strand targetStrand)
	{
		int count = 0;
		while (cpgit.hasNext())
		{
			Cpg cpg = cpgit.next();
			if (cpg.getStrand() == targetStrand) count++; 
		}
		return count;
	}
	
	
}


