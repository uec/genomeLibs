package edu.usc.epigenome.genomeLibs;

import java.util.*;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.*;

public class ReadPos implements Cloneable {
	
	/* Class vars */
	public static final int UNKNOWN = -1; 
	
	
	/* Obj vars */
	protected Symbol sym = DNATools.n();
	protected StrandedFeature.Strand strand = StrandedFeature.UNKNOWN;
	protected static final int STRINGBUFLEN = 10000;
	protected static StringBuffer STRINGBUF = new StringBuffer(STRINGBUFLEN);
	
	
	/* Constructors */
	public ReadPos()
	{
	}

	public ReadPos(Symbol inSym, boolean inForwardStrand)
	{
		this.sym = inSym;
		this.strand = (inForwardStrand) ? StrandedFeature.POSITIVE : StrandedFeature.NEGATIVE;
	}

	/* Getters/Setters */

	/**
	 * @return the symbol relative to the orientation of the genome
	 */
	public Symbol getSym() {
		return sym;
	}
	
	/**
	 * @return the symbol relative the orientation of the read
	 */
	public Symbol getSymReaddir()
	{
		
		Symbol out = sym;
		
		if (this.getStrand() == StrandedFeature.NEGATIVE)
		{
			try 
			{
				out = DNATools.complement(sym);
			}
			catch (IllegalSymbolException e)
			{
				System.err.println("Can not complement symbol " + sym.getName());
				e.printStackTrace();
				System.exit(0);
			}
		}
		
		return out;
	}
	
	public char getSymToken()
	{
		Symbol sym = this.getSym();
		char c = '0';
		try
		{
			c = DNATools.dnaToken(sym);
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(0);
		}
		
		return c;
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
	 * @param strand the strand to set
	 */
	public void setStrand(StrandedFeature.Strand strand) {
		this.strand = strand;
	}


	/**
	 * @return the readPos
	 */
	public int getCycle() {
		return UNKNOWN;
	}


	/**
	 * @return the qual
	 */
	public int getQual() {
		return UNKNOWN;
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
			ReadPos anotherReadPos = (ReadPos)obj;
			String thisKey = this.equalsKey();
			String thatKey = anotherReadPos.equalsKey();
			//System.err.println("Comparing key " + thisKey + " to " + thatKey);
			return thisKey.equals(thatKey);
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
		
	//TODO Not threadsafe because of static buffer STRINGBUF
	public String commaSeparatedLine()
	{
		char delim = ',';		
		
		if (true)
		{
//			System.err.println("New method");
			// New, this takes about 35% of total execution time

			STRINGBUF.delete(0, STRINGBUFLEN);
			STRINGBUF.append(this.getSymReaddir().getName());
			STRINGBUF.append(delim);
			STRINGBUF.append(this.getStrand());
			STRINGBUF.append(delim);
			STRINGBUF.append(this.getCycle());
			STRINGBUF.append(delim);
			STRINGBUF.append(this.getQual());
			return STRINGBUF.toString();
		}
		else
		{
//			System.err.println("Old method");
			// old.  67% of total execution time was spend in string concatenation
			String key = "";
			key += this.getSymReaddir().getName();
			key += delim;
			key += this.getStrand();
			key += delim;
			key += this.getCycle();
			key += delim;
			key += this.getQual();
			return key;
		}
	}
	

}
