package edu.usc.epigenome.genomeLibs.MethylDb;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.seq.StrandedFeature;

public class Cytosine implements Comparable, Cloneable {

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
	
	public short preBaseTotalReads = 0;
	public short preBaseGreads = 0;
	
	// This is actually for non-cpgs
	protected char nextBaseRefUpperCase = '0';
	protected char preBaseRefUpperCase = '0';
	
	protected double methyDens = Double.NaN;
	
	// This is for weighted averages
	protected double gchWeight = Double.NaN;
	protected double gcgWeight = Double.NaN;
	protected double hcgWeight = Double.NaN;



	/**
	 * 
	 */
	public Cytosine() {
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
	public Cytosine(int chromPos, boolean negStrand) {
		super();
		this.chromPos = chromPos;
		this.negStrand = negStrand;
	}
	
	public Cytosine(int chromPos, boolean negStrand, short totalReads, short cReads,
			short cReadsNonconversionFilt, short tReads, short agReads,
			short totalReadsOpposite, short aReadsOpposite, String preBaseRefUpperCase, String nextBaseRefUpperCase, double fracMeth, int gchWeight,  int gcgWeight, int hcgWeight ) {
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
		
		setPreBaseRef(preBaseRefUpperCase);
		setNextBaseRef(nextBaseRefUpperCase);
		
		this.methyDens = fracMeth;
		
		// CpG weight has some issues.  For instance at large gaps you get huge values.  Remove these
		if (gchWeight < 0)
		{
			System.err.printf("Got negative gchWeight: pos=%d, weight=%d\n",chromPos,gchWeight);
			gchWeight = 100;
		}
		else if (gchWeight > 5000)
		{
			System.err.printf("Got extra-large gchWeight: pos=%d, weight=%d\n",chromPos,gchWeight);
			gchWeight = 100;
		}
		
		if (gcgWeight < 0)
		{
			System.err.printf("Got negative gcgWeight: pos=%d, weight=%d\n",chromPos,gcgWeight);
			gcgWeight = 100;
		}
		else if (gcgWeight > 5000)
		{
			System.err.printf("Got extra-large gcgWeight: pos=%d, weight=%d\n",chromPos,gcgWeight);
			gcgWeight = 100;
		}
		
		if (hcgWeight < 0)
		{
			System.err.printf("Got negative hcgWeight: pos=%d, weight=%d\n",chromPos,hcgWeight);
			hcgWeight = 100;
		}
		else if (hcgWeight > 5000)
		{
			System.err.printf("Got extra-large hcgWeight: pos=%d, weight=%d\n",chromPos,hcgWeight);
			hcgWeight = 100;
		}
		
		// It's actually divided by two since we're looking at both strands.
		this.gchWeight = (double)gchWeight/2.0;
		this.gcgWeight = (double)gcgWeight/2.0;
		this.hcgWeight = (double)hcgWeight/2.0;
	}
	
	/*** Overridden Comparable methods
	 */

	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Object o) {

		Cytosine other = (Cytosine)o;
		
		return (this.chromPos - other.chromPos);
	}


	/* (non-Javadoc)
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object obj) {
		return 	(this.compareTo(obj) == 0);
	}
	
	public static Cytosine compositeCytocine(List<Cytosine> inCytocines) throws CloneNotSupportedException
	{
		Cytosine out = null;
		
		for (Cytosine inCytocine : inCytocines)
		{
			if (out == null)
			{
				// First one
				out = (Cytosine)inCytocine.clone();
			}
			else
			{
				out.totalReads += inCytocine.totalReads;
				out.cReads += inCytocine.cReads;
				out.cReadsNonconversionFilt += inCytocine.cReadsNonconversionFilt;
				out.tReads += inCytocine.tReads;
				out.agReads += inCytocine.agReads;
				out.totalReadsOpposite += inCytocine.totalReadsOpposite;
				out.aReadsOpposite += inCytocine.aReadsOpposite;
			}
		
		}
		
		return out;
	}
	
	public double fracNextBaseG()
	{
		return (double)this.nextBaseGreads / ((double)this.nextBaseTotalReads);
	}
	
	public double fracPreBaseG()
	{
		return (double)this.preBaseGreads / ((double)this.preBaseTotalReads);
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

	public double fracMeth(){
		return methyDens;
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

	
	public static void outputCytocinesToFile(PrintWriter pw, Map<Integer,Cytosine> cytocineMap, String prefix, String sampleName, String chr)
	throws IOException
	{
		//System.err.println("About to write " + cpgMap.size() + " Cytosines to file");
		Iterator<Cytosine> cytocineIt = cytocineMap.values().iterator();
		 while (cytocineIt.hasNext())
		{
			Cytosine cytocine = cytocineIt.next();
			//if ((cytocine.getNextBaseRef() != 'G'))// && (cpg.getNextBaseRef() != '0'))
			//{
				
				//	continue CPG;
				
			//}
			//if ((cytocine.getLastBaseRef() == 'G') && (cytocine.getNextBaseRef() != 'G')){
				
			//}
			//else if((cytocine.getLastBaseRef() == 'G') && (cytocine.getNextBaseRef() == 'G')){
				
			//}
			//else if((cytocine.getLastBaseRef() != 'G') && (cytocine.getNextBaseRef() == 'G')){
				
			//}
			//else{
				
			//}
			
			String line = cytocine.toString();
			pw.println(line);
		}
	}
	
	public static PrintWriter outputChromToFile(Map<Integer,Cytosine> cytocineMap, String prefix, String sampleName, String chr)
	throws IOException
	{
		
		String fn = prefix + sampleName + "_" + chr + ".txt";
		PrintWriter writer = new PrintWriter(new File(fn));
		
		outputCytocinesToFile(writer, cytocineMap, prefix, sampleName, chr);
		
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

		return String.format("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c\t%d\t%d\t%c\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f", 
				chromPos,
				(negStrand) ? '-' : '+',
				totalReads,
				cReads,
				cReadsNonconversionFilt,
				tReads,
				agReads,
				totalReadsOpposite,
				aReadsOpposite,
				preBaseGreads,
				preBaseTotalReads,
				preBaseRefUpperCase,
				nextBaseGreads,
				nextBaseTotalReads,
				nextBaseRefUpperCase,
				100*this.fracMeth(true),
				100*this.fracPreBaseG(),
				100*this.fracNextBaseG(),
				this.gchWeight,
				this.gcgWeight,
				this.hcgWeight
				);
				
	}
	

	@Override
	protected Object clone() throws CloneNotSupportedException {
		Cytosine out = (Cytosine)super.clone();
		
		System.err.printf("Cloning cytocine, old tReads=%d, new tReads=%d\n",  this.tReads, out.tReads);
		return out;
	}

	public double getGchWeight() {
		return gchWeight;
	}
	
	public double getGcgWeight() {
		return gcgWeight;
	}
	
	public double getHcgWeight() {
		return hcgWeight;
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
	
	public void setPreBaseRef(String preBaseRef){
		this.preBaseRefUpperCase = preBaseRef.charAt(0);
	}
	
	public void setNextBaseRef(String nextBaseRef){
		this.nextBaseRefUpperCase = nextBaseRef.charAt(0);
	}
	
	public char getPreBaseRef() {
		return preBaseRefUpperCase;
	}

	public void setPreBaseRef(char preBaseRef) {
		this.preBaseRefUpperCase = Character.toUpperCase(preBaseRef);
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

	public Cytosine downsample(double downsamplingFactor) {
		// TODO Auto-generated method stub
		
		System.err.println("Cytocine::downsample not yet supported");
		System.exit(1);
		
		return null;
	}
	
	
	
	
}