package edu.usc.epigenome.genomeLibs.MethylDb;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Map;

import org.biojava.bio.seq.StrandedFeature;

public class Cpg implements Comparable {

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
	protected char nextBaseRefUpperCase = '0';
	
	// This is for weighted averages
	protected double cpgWeight = Double.NaN;



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
	
	public Cpg(int chromPos, boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite) {
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

	
	public double fracMeth(boolean useNonconvFilt)
	{
		int cs = (useNonconvFilt) ? this.cReads : (this.cReads + this.cReadsNonconversionFilt); 
		return (double)cs / ((double)cs + (double)this.tReads);
	}
	
	public String variableStepWigLine(boolean useNonconvFilter)
	{
		return String.format("%d\t%f",chromPos, this.fracMeth(useNonconvFilter));
	}

	public static void outputChromToFile(Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr)
	throws IOException
	{
		outputChromToFile(cpgMap, prefix, sampleName, chr, 0, 0.0);
	}
	
	public static void outputCpgsToFile(PrintWriter pw, Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr, int minCphCoverage, double minCphMethFrac)
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
			
			String line = cpg.toString();
			pw.println(line);
		}
	}
	
	public static PrintWriter outputChromToFile(Map<Integer,Cpg> cpgMap, String prefix, String sampleName, String chr, int minCphCoverage, double minCphMethFrac)
	throws IOException
	{
		
		String fn = prefix + sampleName + "_" + chr + ".txt";
		PrintWriter writer = new PrintWriter(new File(fn));
		
		outputCpgsToFile(writer, cpgMap, prefix, sampleName, chr, minCphCoverage, minCphMethFrac);
		
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

		return String.format("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%.2f\t%.2f", 
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
				100*this.fracMeth(true),
				100*this.fracNextBaseG()
				);
	}
	

	public double getCpgWeight() {
		return cpgWeight;
	}

	public void setCpgWeight(double cpgWeight) {
		this.cpgWeight = cpgWeight;
	}

	public char getNextBaseRef() {
		return nextBaseRefUpperCase;
	}

	public void setNextBaseRef(char nextBaseRef) {
		this.nextBaseRefUpperCase = Character.toUpperCase(nextBaseRef);
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
		
		System.err.println("Cpg::downsample not yet supported");
		System.exit(1);
		
		return null;
	}
	
	
	
	
}
