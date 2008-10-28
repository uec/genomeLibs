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
	 * @return the sym
	 */
	public Symbol getSym() {
		return sym;
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
	public int getReadPos() {
		return UNKNOWN;
	}


	/**
	 * @return the qual
	 */
	public int getQual() {
		return UNKNOWN;
	}
	
	
	/* Static util functions */
	
	/**
	 * @param a list of ReadPos objs
	 * @return int[0]=fw_depth , int[1]=rev_depth
	 * 
	 * We only make it non-static so that it can be overridden
	 */

	public int[] getDepth(Collection<ReadPos> posList, ReadPosOptions ro)
	{
		int[] depth = new int[] {0,0};

		Iterator<ReadPos> it = posList.iterator();
		while (it.hasNext())
		{
			ReadPos rp = it.next();
			int index = (rp.getStrand()==StrandedFeature.NEGATIVE) ? 1 : 0;
			depth[index]++;
		}
		
		return depth;
	}

	
	
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
	
	

}
