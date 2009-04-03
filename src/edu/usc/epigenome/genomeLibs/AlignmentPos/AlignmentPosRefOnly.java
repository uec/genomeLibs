package edu.usc.epigenome.genomeLibs.AlignmentPos;


import java.util.Vector;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.*;


public class AlignmentPosRefOnly extends AlignmentPos {

	/*****************
	 *  Constructors
	 */


	
	public AlignmentPosRefOnly(char inRef, String inChr, int inPos) {
		super(inRef, inChr, inPos, new AlignmentPosOptions());
	}

	public AlignmentPosRefOnly(Symbol inRef, String inChr, int inPos) {
		super(inRef, inChr, inPos, new AlignmentPosOptions());
	}

	public AlignmentPosRefOnly(AlignmentPos ap)
	{
		super(ap);
	}
	
	/**
	 * 
	 * Getters
	 * 
	 */
	
	public int[] getDepth()
	{
		int[] out = new int[] {0,0};
		
		int idx = (this.getStrand()==StrandedFeature.NEGATIVE) ? 1 : 0;
		out[idx]++;
		
		return out;
	}

	
	public  AlignmentPosRefOnly clone(boolean flipStrand)
	{
		AlignmentPosRefOnly ap = new AlignmentPosRefOnly(this.getRefFlipped(), this.chr, this.pos);
		StrandedFeature.Strand strand = (flipStrand) ? this.getStrand().flip() : this.getStrand();
		ap.setStrand(strand);

		return ap;
	}

	

	/*****************
	 * 
	 * Setters
	 * 
	 */
	
	
	public void removeRevStrandReads()
	{
	}
	
	public void reset()
	{
	}
	


}
