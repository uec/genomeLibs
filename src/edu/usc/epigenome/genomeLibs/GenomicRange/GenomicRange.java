package edu.usc.epigenome.genomeLibs.GenomicRange;

import org.biojava.bio.seq.StrandedFeature;


public class GenomicRange implements Cloneable, Comparable<GenomicRange> {
	
    
	//TODO Use org.biojavax.bio.seq.SimplePosition for position info
	
	/* Class vars */
	static final String NO_CHROM = "NO_CHROM";
	static final double NO_SCORE = Double.NEGATIVE_INFINITY;
	
	
	/* Obj vars */
	protected String chrom = NO_CHROM;
	protected int start = 0;
	protected int end = 0;
	protected StrandedFeature.Strand strand = StrandedFeature.UNKNOWN;
	protected double score = NO_SCORE;
	
	/* Constructors */
	public GenomicRange(String inChrom, int inStart, int inEnd)
	{
		chrom = inChrom;
		start = inStart;
		end = inEnd;
	}

	public GenomicRange(String inChrom, int inStart, int inEnd, StrandedFeature.Strand inStrand)
	{
		chrom = inChrom;
		start = inStart;
		end = inEnd;
		strand = inStrand;
	}

	
	/**
	 * @return the strand
	 */
	public StrandedFeature.Strand getStrand() {
		return strand;
	}




	/**
	 * @return the chrom
	 */
	public String getChrom() {
		return chrom;
	}



	/**
	 * @param chrom the chrom to set
	 */
	public void setChrom(String chrom) {
		this.chrom = chrom;
	}



	/**
	 * @return the start
	 */
	public int getStart() {
		return start;
	}



	/**
	 * @param strand the strand to set
	 */
	public void setStrand(StrandedFeature.Strand strand) {
		this.strand = strand;
	}


	/**
	 * @param start the start to set
	 */
	public void setStart(int start) {
		this.start = start;
	}



	/**
	 * @return the end
	 */
	public int getEnd() {
		return end;
	}



	/**
	 * @param end the end to set
	 */
	public void setEnd(int end) {
		this.end = end;
	}


	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}



	/**
	 * @param score the score to set
	 */
	public void setScore(double score) {
		this.score = score;
	}



	/** Static util functions **/
	public static GenomicRange generateWindowFromCoord(String inChrom, int inCoord, int inWindSize)
	{
			return generateWindowFromCoord(inChrom, inCoord, inWindSize, StrandedFeature.UNKNOWN);
	}	
	
	public static GenomicRange generateWindowFromCoord(String inChrom, int inCoord, int inWindSize, StrandedFeature.Strand strand)
	{
		int s = ((int)Math.floor(inCoord/inWindSize))*inWindSize;
		int e = s + inWindSize - 1;
		return new GenomicRange(inChrom, s, e, strand);
	}
	

	/** Cloneable **/
	
	public GenomicRange clone()
	{
		GenomicRange newInst = null;
		try
		{
			newInst = (GenomicRange)super.clone();
		}
		catch (CloneNotSupportedException e)
		{
			// Should never get here
			e.printStackTrace();
			System.exit(0);
		}
		
		return newInst;
	}

	/** 
	 * If symbol and strand are equivalent, the two are equal
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj) return true;
		
		if (obj instanceof GenomicRange)
		{
			return (this.compareTo((GenomicRange)obj) == 0);
		}
		else
		{	
			return false;
		}
	}

	/**
	 * If symbol and strand are equivalent, the two are equal
	 */
	@Override
	public int hashCode() {
		return this.equalsKey().hashCode();
	}
	
	protected String equalsKey()
	{
		return commaSeparatedLine();
	}
		
	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(GenomicRange o) {
		
		int strCmp = this.getChrom().compareTo(o.getChrom());
		if (strCmp != 0) return strCmp;

		// Do the end first (this only matters if one completely 
		// subsumes the other, in which case the internal one will
		// be returned first.
		if (this.getEnd() > o.getEnd()) return 1; 
		if (this.getEnd() < o.getEnd()) return -1;
		
		if (this.getStart() > o.getStart()) return 1; 
		if (this.getStart() < o.getStart()) return -1;
		
		if (strandToInt(this.getStrand()) > strandToInt(o.getStrand())) return 1; 
		if (strandToInt(this.getStrand()) < strandToInt(o.getStrand())) return -1; 

		return 0;
	}

	
	public static int strandToInt(StrandedFeature.Strand strand)
	{
		int out = 0;
		
		if (strand == StrandedFeature.NEGATIVE)
		{
			out = -1;
		}
		else if (strand == StrandedFeature.POSITIVE)
		{
			out = 1;
		}
		
		return out;
	}
	
	public String commaSeparatedLine()
	{
		return String.format("%s,%d,%d,%d",this.getChrom(), strandToInt(this.getStrand()), this.getStart(), this.getEnd());
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return this.commaSeparatedLine();
	}
	

	
	
}
