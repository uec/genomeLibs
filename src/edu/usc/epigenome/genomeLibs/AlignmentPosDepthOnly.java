package edu.usc.epigenome.genomeLibs;


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

	
	public  AlignmentPosDepthOnly clone(boolean flip_strand)
	{
		AlignmentPosDepthOnly ap = new AlignmentPosDepthOnly(this.getRefFlipped(), this.chr, this.pos, this.apOptions);
		ap.setDepth(this.getDepth());
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

	
	public void reset()
	{
		this.depths = new int[] {0,0};
	}
	


}
