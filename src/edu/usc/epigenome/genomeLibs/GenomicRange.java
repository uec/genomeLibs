package edu.usc.epigenome.genomeLibs;


public class GenomicRange implements Cloneable, Comparable<GenomicRange> {
	

	/* Obj vars */
	protected String chrom = null;
	protected int start = 0;
	protected int end = 0;
	
	/* Constructors */
	public GenomicRange(String inChrom, int inStart, int inEnd)
	{
		chrom = inChrom;
		start = inStart;
		end = inEnd;
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


	/** Static util functions **/
	public static GenomicRange generateWindowFromCoord(String inChrom, int inCoord, int inWindSize)
	{
		int s = ((int)Math.floor(inCoord/inWindSize))*inWindSize;
		int e = s + inWindSize - 1;
		return new GenomicRange(inChrom, s, e);
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

		return 0;
	}

	public String commaSeparatedLine()
	{
		return String.format("%s,%d,%d",this.getChrom(), this.getStart(), this.getEnd());
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return this.commaSeparatedLine();
	}
	

	
	
}
