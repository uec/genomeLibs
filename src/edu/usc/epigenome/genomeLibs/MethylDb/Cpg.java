package edu.usc.epigenome.genomeLibs.MethylDb;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.DistributionFactory;
import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.MatUtils;

public class Cpg implements Comparable, Cloneable {

	private static final double DEFAULT_FAIL_MIN_FRAC_A = 0.2;
	private static final int DEFAULT_FAIL_MIN_ABS_A = 2;
	
	// Object vars
	public int chromPos = 0;
	public boolean negStrand = false;
	public short totalReads = 0;
	public short cReads = 0;
	public short cReadsNonconversionFilt = 0;
	public short tReads = 0;
	public short agReads = 0;
	public short totalReadsOpposite = 0;
	public short aReadsOpposite = 0;
	
	public short nextBaseTotalReads = 0;
	public short nextBaseGreads = 0;
	
	// This is actually for non-cpgs
	protected char prevBaseRefUpperCase = '0';
	protected char nextBaseRefUpperCase = '0';
	
	// This is for weighted averages
	protected double cpgWeight = Double.NaN;
	
	// Read info
	protected Map<Integer,CpgRead> readMap = new TreeMap<Integer,CpgRead>();



	/**
	 * 
	 */
	public Cpg() {
		super();
	}

	/**
	 * @param chromPos
	 * @param negStrand
	 * @param totalReads
	 * @param cReads
	 * @param cReadsNonconversionFilt
	 * @param tReads
	 * @param agReads
	 * @param totalReadsOpposite
	 * @param aReadsOpposite
	 */
	public Cpg(int chromPos, boolean negStrand) {
		super();
		this.chromPos = chromPos;
		this.negStrand = negStrand;
	}
	
//	public Cpg(int chromPos, boolean negStrand, short totalReads, short cReads,
//			short cReadsNonconversionFilt, short tReads, short agReads,
//			short totalReadsOpposite, short aReadsOpposite, int cpgWeight) {
//		super();
//		this.chromPos = chromPos;
//		this.negStrand = negStrand;
//		this.totalReads = totalReads;
//		this.cReads = cReads;
//		this.cReadsNonconversionFilt = cReadsNonconversionFilt;
//		this.tReads = tReads;
//		this.agReads = agReads;
//		this.totalReadsOpposite = totalReadsOpposite;
//		this.aReadsOpposite = aReadsOpposite;
//		
//		// CpG weight has some issues.  For instance at large gaps you get huge values.  Remove these
//		if (cpgWeight < 0)
//		{
//			System.err.printf("Got negative cpgWeight: pos=%d, weight=%d\n",chromPos,cpgWeight);
//			cpgWeight = 100;
//		}
//		else if (cpgWeight > 5000)
//		{
//			System.err.printf("Got extra-large cpgWeight: pos=%d, weight=%d\n",chromPos,cpgWeight);
//			cpgWeight = 100;
//		}
//		
//		// It's actually divided by two since we're looking at both strands.
//		this.cpgWeight = (double)cpgWeight/2.0;
//	}
	
	public Cpg(int chromPos, boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite, int cpgWeight,
			short nextBaseGreads, short nextBaseTotalReads, char nextBaseRefUpperCase)
	 {
		super();
		this.chromPos = chromPos;
		this.negStrand = negStrand;
		this.totalReads = totalReads;
		this.cReads = cReads;
		this.cReadsNonconversionFilt = cReadsNonconversionFilt;
		this.tReads = tReads;
		this.agReads = agReads;
		this.totalReadsOpposite = totalReadsOpposite;
		this.aReadsOpposite = aReadsOpposite;
		
		this.nextBaseGreads = nextBaseGreads;
		this.nextBaseTotalReads = nextBaseTotalReads;
		this.nextBaseRefUpperCase = nextBaseRefUpperCase;
		
		// CpG weight has some issues.  For instance at large gaps you get huge values.  Remove these
		if (cpgWeight < 0)
		{
			System.err.printf("Got negative cpgWeight: pos=%d, weight=%d\n",chromPos,cpgWeight);
			cpgWeight = 100;
		}
		else if (cpgWeight > 5000)
		{
			System.err.printf("Got extra-large cpgWeight: pos=%d, weight=%d\n",chromPos,cpgWeight);
			cpgWeight = 100;
		}
		
		// It's actually divided by two since we're looking at both strands.
		this.cpgWeight = (double)cpgWeight/2.0;
	}
	
	/*** Read information ***/
	public void addRead(CpgRead read)
	{
		readMap.put(read.readId, read);
		
		this.agReads += read.agRead;
		this.tReads += read.tRead;
		this.cReads += read.cRead;
		this.totalReads += (read.agRead + read.tRead + read.cRead);
		this.cReadsNonconversionFilt += read.cReadNonconversionFilt;
		this.nextBaseGreads += read.nextBaseGread;
	}
	
	public Map<Integer,CpgRead> getReads()
	{
		return readMap;
	}
	
	
	/*** Overridden Comparable methods
	 */

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Object o) {

		Cpg other = (Cpg)o;
		
		return (this.chromPos - other.chromPos);
	}


	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		return 	(this.compareTo(obj) == 0);
	}
	
	public static Cpg compositeCpg(List<Cpg> inCpgs) throws CloneNotSupportedException
	{
		Cpg out = null;
		
		for (Cpg inCpg : inCpgs)
		{
			if (out == null)
			{
				// First one
				out = (Cpg)inCpg.clone();
			}
			else
			{
				out.totalReads += inCpg.totalReads;
				out.cReads += inCpg.cReads;
				out.cReadsNonconversionFilt += inCpg.cReadsNonconversionFilt;
				out.tReads += inCpg.tReads;
				out.agReads += inCpg.agReads;
				out.totalReadsOpposite += inCpg.totalReadsOpposite;
				out.aReadsOpposite += inCpg.aReadsOpposite;
			}
		
		}
		
		return out;
	}
	
	public double fracNextBaseG()
	{
		return (double)this.nextBaseGreads / ((double)this.nextBaseTotalReads);
	}

	public double fracOppositeA()
	{
		return (double)this.aReadsOpposite / ((double)this.aReadsOpposite + (double)this.totalReadsOpposite);
	}
	
	public boolean passesOppositeAFilterDefault()
	{
		return passesOppositeAFilter(DEFAULT_FAIL_MIN_FRAC_A, DEFAULT_FAIL_MIN_ABS_A);
	}
	
	public boolean passesOppositeAFilter(double failMinFracA, int failMinAbsA)
	{
		double fracOppA = this.fracOppositeA();
		boolean failFrac = (!Double.isNaN(fracOppA)) && (fracOppA>=failMinFracA);
		boolean failAbs = (this.aReadsOpposite>=failMinAbsA);
		return !(failFrac && failAbs);
	}

	public double fracReadsCorT()
	{
		return (double)totalReadsCorT(false) / ((double)totalReadsCorT(false)+(double)agReads);
	}
	
	public double fracMeth(boolean useNonconvFilt)
	{
		int cs = (useNonconvFilt) ? this.cReads : (this.cReads + this.cReadsNonconversionFilt); 
		return (double)cs / ((double)cs + (double)this.tReads);
	}
	
	public String variableStepWigLine(boolean useNonconvFilter)
	{
		return String.format("%d\t%f",chromPos, this.fracMeth(useNonconvFilter));
	}

	public static void outputCpgsToFile(PrintWriter pw, Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr, int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		outputCpgsToFile(false, pw, cpgMap, prefix, sampleName, chr, minCphCoverage, minCphMethFrac);
	}
	
	public static void outputCpgsToFile(boolean readFormat, PrintWriter pw, Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr, int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		//System.err.println("About to write " + cpgMap.size() + " Cytosines to file");
		Iterator<Cpg> cpgIt = cpgMap.values().iterator();
		CPG: while (cpgIt.hasNext())
		{
			Cpg cpg = cpgIt.next();
			if ((cpg.getNextBaseRef() != 'G'))// && (cpg.getNextBaseRef() != '0'))
			{
				if ((cpg.totalReads < minCphCoverage) || (cpg.fracMeth(true)<minCphMethFrac))
				{
					continue CPG;
				}
			}
			
			String line = (readFormat) ? cpg.toReadFormatString() : cpg.toString();
			pw.print(line);
			if (!readFormat) pw.println();
		}
	}
	


	public static PrintWriter outputChromToFile(Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr, int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		return outputChromToFile(false, cpgMap, prefix, sampleName, chr,  minCphCoverage, minCphMethFrac);
	}
		
	public static void outputChromToFile(Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr)
	throws IOException
	{
		outputChromToFile(cpgMap, prefix, sampleName, chr, 0, 0.0);
	}
	
	public static PrintWriter outputChromToFile(boolean readFormat, Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr, 
			int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		
		String fn = prefix + sampleName + "_" + chr + ".txt";
		PrintWriter writer = new PrintWriter(new File(fn));
		
		outputCpgsToFile(readFormat, writer, cpgMap, prefix, sampleName, chr, minCphCoverage, minCphMethFrac);
		
		//writer.close();
		return writer;
	}

	@Override
	public String toString() {

		if (true)
		{
			return this.toStringExpanded();
		}
		else
		{
		return String.format("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d", 
				chromPos,
				(negStrand) ? '-' : '+',
				totalReads,
				cReads,
				cReadsNonconversionFilt,
				tReads,
				agReads,
				totalReadsOpposite,
				aReadsOpposite
				);
		}
	}
	
	public String toStringExpanded() 
	{

		return String.format("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%c\t%.2f\t%.2f\t%.2f", 
				chromPos,
				(negStrand) ? '-' : '+',
				totalReads,
				cReads,
				cReadsNonconversionFilt,
				tReads,
				agReads,
				totalReadsOpposite,
				aReadsOpposite,
				nextBaseGreads,
				nextBaseTotalReads,
				nextBaseRefUpperCase,
				prevBaseRefUpperCase,
				100*this.fracMeth(true),
				100*this.fracNextBaseG(),
				this.cpgWeight
				);
	}
	
	public String toReadFormatString() {
		
		StringBuffer sb = new StringBuffer(100*this.readMap.size());
		Iterator<CpgRead> readIt = this.readMap.values().iterator();
		while (readIt.hasNext())
		{
			CpgRead read = readIt.next();
			

			sb.append(String.format("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%.2f\n", 
					chromPos,
					((negStrand) ? '-' : '+'),
					read.readId,
					
					read.cRead,
					read.cReadNonconversionFilt,
					read.tRead,
					read.agRead,
					
					totalReadsOpposite,
					aReadsOpposite,
					
					read.nextBaseGread,
					read.nextBaseUpperCase,
					
					this.cpgWeight
			));
		}

//		  `chromPos` INT UNSIGNED NOT NULL,
//		  `strand` enum('+','-') NOT NULL,
//		  `readId` INT UNSIGNED	NOT NULL,
//
//		  `cRead` SMALLINT UNSIGNED NOT	NULL,
//		  `cReadNonconversionFilt` SMALLINT UNSIGNED NOT NULL,
//		  `tRead` SMALLINT UNSIGNED NOT NULL,
//		  `agRead` SMALLINT UNSIGNED NOT NULL,
//
//		  `totalReadsOpposite` SMALLINT UNSIGNED NOT NULL,
//		  `aReadsOpposite` SMALLINT UNSIGNED NOT NULL,
//
//		  `nextBaseGread` SMALLINT UNSIGNED NULL default '0',
//		  `cpgWeight` SMALLINT UNSIGNED NOT NULL,

		return sb.toString();
	}
	
	
	@Override
	protected Object clone() throws CloneNotSupportedException {
		Cpg out = (Cpg)super.clone();
		
		//System.err.printf("Cloning cpg, old tReads=%d, new tReads=%d\n",  this.tReads, out.tReads);
		return out;
	}

	public double getCpgWeight() {
		return cpgWeight;
	}

//	public void setCpgWeight(double cpgWeight) {
//		this.cpgWeight = cpgWeight;
//	}

	public char getNextBaseRef() {
		return nextBaseRefUpperCase;
	}

	public void setNextBaseRef(char nextBaseRef) {
		this.nextBaseRefUpperCase = Character.toUpperCase(nextBaseRef);
	}

	/**
	 * @return the prevBaseRefUpperCase
	 */
	public char getPrevBaseRef() {
		return prevBaseRefUpperCase;
	}

	/**
	 * @param prevBaseRefUpperCase the prevBaseRefUpperCase to set
	 */
	public void setPrevBaseRef(char prevBaseRefUpperCase) {
		this.prevBaseRefUpperCase = prevBaseRefUpperCase;
	}
	
	public StrandedFeature.Strand getStrand()
	{
		return (this.negStrand) ? StrandedFeature.NEGATIVE : StrandedFeature.POSITIVE;
	}
	
	public String getStrandStr()
	{
		return (this.negStrand) ? "-" : "+";
	}
	
	public short totalReadsC(boolean useNonconvFilt)
	{
		int cs = (useNonconvFilt) ? this.cReads : (this.cReads + this.cReadsNonconversionFilt); 
		return (short)cs;
	}
	
	public short totalReadsCorT(boolean useNonconvFilt)
	{
		short cs = totalReadsC(useNonconvFilt);
		return (short)(cs + this.tReads);
	}

	public Cpg downsample(double downsamplingFactor) {
		// TODO Auto-generated method stub
		
		
		Cpg out = null;
		try
		{
			out = (Cpg)this.clone();
			
			// For positive strand, we have all the counts so we can sample every flavor separately
			out.cReads = (short)MatUtils.randomBinomialGeneratedCount(this.cReads, downsamplingFactor);
			out.cReadsNonconversionFilt = (short)MatUtils.randomBinomialGeneratedCount(this.cReadsNonconversionFilt, downsamplingFactor);
			out.tReads = (short)MatUtils.randomBinomialGeneratedCount(this.tReads, downsamplingFactor);
			out.agReads = (short)MatUtils.randomBinomialGeneratedCount(this.agReads, downsamplingFactor);
			out.totalReads = (short)(out.cReads + out.cReadsNonconversionFilt + out.tReads + out.agReads);
			
			// For reverse strand and next base, we only have counts for one variant, so we have to resample the total
			// and then sample the single variant based on the observed fraction
			out.totalReadsOpposite = (short)MatUtils.randomBinomialGeneratedCount(this.totalReadsOpposite, downsamplingFactor);
			out.aReadsOpposite = (short)MatUtils.randomBinomialGeneratedCount(out.totalReadsOpposite, this.fracOppositeA());
			out.nextBaseTotalReads = (short)MatUtils.randomBinomialGeneratedCount(this.nextBaseTotalReads, downsamplingFactor);
			out.nextBaseGreads = (short)MatUtils.randomBinomialGeneratedCount(out.nextBaseTotalReads, this.fracNextBaseG());
		}
		catch (Exception e)
		{
			System.err.printf("Problem downsampling: %s\n",e.toString());
			e.printStackTrace();
			System.exit(1);
		}

		// - - - REMOVE --- REMOVE --- REMOVE !!!
//		System.out.printf("%d,%.2f,%.2f,%.2f,%.2f\n", this.totalReads, this.fracMeth(true), out.fracMeth(true), this.fracMeth(false), out.fracMeth(false));

		
//		System.err.printf("Before: %s\n",this.toStringExpanded());
//		System.err.printf("After: %s\n\n",out.toStringExpanded());
			
		
		
//		System.err.println("Cpg::downsample not yet supported");
//		System.exit(1);
		
		return out;
	}
	
	
	
	
}
