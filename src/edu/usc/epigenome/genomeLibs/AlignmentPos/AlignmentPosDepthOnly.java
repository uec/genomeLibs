package edu.usc.epigenome.genomeLibs.AlignmentPos;


import java.util.Vector;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.*;


public class AlignmentPosDepthOnly extends AlignmentPos {

	/* Obj vars */
	protected int[] depths = new int[] {0,0};


	/*****************
	 *  Constructors
	 */


	
	public AlignmentPosDepthOnly(char inRef, String inChr, int inPos, AlignmentPosOptions inApOptions) {
		super(inRef, inChr, inPos,inApOptions);
	}


	public AlignmentPosDepthOnly(Symbol inRef, String inChr, int inPos, AlignmentPosOptions inApOptions) {
		super(inRef, inChr, inPos,inApOptions);
	}

	public AlignmentPosDepthOnly(AlignmentPos ap)
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
		return depths;
	}

	
	public  AlignmentPosDepthOnly clone(boolean flipStrand)
	{
		AlignmentPosDepthOnly ap = new AlignmentPosDepthOnly(this.getRefFlipped(), this.chr, this.pos, this.apOptions);
		ap.depths[0] = this.depths[1];
		ap.depths[1] = this.depths[0];
		
		StrandedFeature.Strand strand = (flipStrand) ? this.getStrand().flip() : this.getStrand();
		ap.setStrand(strand);

		return ap;
	}

	

	/*****************
	 * 
	 * Setters
	 * 
	 */
	
	public void setDepth(int[] inDepth)
	{
		this.depths[0] = inDepth[0];
		this.depths[1] = inDepth[1];
	}
	
	public void incrementDepth(int[] inDepth)
	{
		this.depths[0] += inDepth[0];
		this.depths[1] += inDepth[1];
	}
	
	public void removeRevStrandReads()
	{
		this.depths[1] = 0;
	}
	
	public void reset()
	{
		this.depths = new int[] {0,0};
	}
	


}
