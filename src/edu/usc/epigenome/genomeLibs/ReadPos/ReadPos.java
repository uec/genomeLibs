package edu.usc.epigenome.genomeLibs.ReadPos;


import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.*;

import edu.usc.epigenome.genomeLibs.BiojavaUtils;

/**
 * Counted objects, contains symbol, strand info
 * @author benb
 *
 */
public class ReadPos implements Cloneable, Comparable<ReadPos> {
	
	/* Class vars */
	public static final int UNKNOWN_CYCLE = -1; 
	public static final int UNKNOWN_QUAL = -1; 
	
	
	/* Obj vars */
	protected Symbol sym = DNATools.n();
	protected StrandedFeature.Strand strand = StrandedFeature.UNKNOWN;
	
	/* Constructors */
	public ReadPos()
	{
	}

	public ReadPos(Symbol inSym, StrandedFeature.Strand inStrand)
	{
		this.sym = inSym;
		this.strand = inStrand;
		//this.strand = (inForwardStrand) ? StrandedFeature.POSITIVE : StrandedFeature.NEGATIVE;
	}

	/* Getters/Setters */

	/**
	 * @return the symbol relative to the orientation of the genome
	 */
	public Symbol getSym() {
		return sym;
	}
	
//	/**
//	 * @return the symbol relative the orientation of the read
//	 */
//	public Symbol getSymReaddir()
//	{
//		
//		Symbol out = sym;
//		
//		if (this.getStrand() == StrandedFeature.NEGATIVE)
//		{
//			try 
//			{
//				out = DNATools.complement(sym);
//			}
//			catch (IllegalSymbolException e)
//			{
//				System.err.println("Can not complement symbol " + sym.getName());
//				e.printStackTrace();
//				System.exit(0);
//			}
//		}
//		
//		return out;
//	}
	
	public char getSymToken()
	{
		Symbol sym = this.getSym();
		return BiojavaUtils.dnaTokenNoException(sym);
	}
	



	/**
	 * @param sym the sym to set
	 */
	public void setSym(Symbol sym) {
		this.sym = sym;
	}


	/**
	 * @return the strand
	 */
	public StrandedFeature.Strand getStrand() {
		return strand;
	}


	/**
	 * @return the strand
	 */
	public char getStrandChar() {
		return strand.getToken();
	}

	/**
	 * @param strand the strand to set
	 */
	public void setStrand(StrandedFeature.Strand strand) {
		this.strand = strand;
	}


	/**
	 * @return the readPos
	 */
	public int getCycle() {
		return UNKNOWN_CYCLE;
	}


	/**
	 * @return the qual
	 */
	public int getQual() {
		return UNKNOWN_QUAL;
	}
	
	
	public ReadPos reverseComplement()
	{
		ReadPos out = this.clone();
		try
		{
			out.setSym(DNATools.complement(this.getSym()));
		}
		catch (IllegalSymbolException e)
		{
			// Should never get here.
			e.printStackTrace();
			System.exit(0);
		}
		
		out.setStrand(this.getStrand().flip());
		return out;
	}
	
	/* Static util functions */
	
	
	
	/** Cloneable **/
	
	public ReadPos clone()
	{
		ReadPos newInst = null;
		try
		{
			newInst = (ReadPos)super.clone();
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
		
		if (obj instanceof ReadPos)
		{
			return (this.compareTo((ReadPos)obj) == 0);
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
	public int compareTo(ReadPos o) {
		
		if (this.getSymToken() > o.getSymToken()) return 1; 
		if (this.getSymToken() < o.getSymToken()) return -1;
		
		if (this.getStrandChar() > o.getStrandChar()) return 1; 
		if (this.getStrandChar() < o.getStrandChar()) return -1;
		
		if (this.getCycle() > o.getCycle()) return 1; 
		if (this.getCycle() < o.getCycle()) return -1;
		
		if (this.getQual() > o.getQual()) return 1; 
		if (this.getQual() < o.getQual()) return -1;
		
		return 0;
	}

	public String commaSeparatedLine()
	{
		return String.format("%c,%c,%d,%d",getSymToken(), getStrandChar(),
				getCycle(), getQual());
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return this.commaSeparatedLine();
	}
	

	
	
}
