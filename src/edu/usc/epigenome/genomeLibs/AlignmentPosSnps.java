package edu.usc.epigenome.genomeLibs;

import java.util.*;

import org.biojava.bio.symbol.*;

public class AlignmentPosSnps extends AlignmentPos {

	/* Obj vars */
	protected Vector<ReadPos> readPosList = new Vector<ReadPos>();
	


	/*****************
	 *  Constructors
	 */


	
	public AlignmentPosSnps(char inRef, String inChr, int inPos, AlignmentPosOptions inApOptions) {
		super(inRef, inChr, inPos,inApOptions);
	}


	public AlignmentPosSnps(Symbol inRef, String inChr, int inPos, AlignmentPosOptions inApOptions) {
		super(inRef, inChr, inPos,inApOptions);
	}

	public AlignmentPosSnps(AlignmentPos ap)
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
		int [] out = new int[] {0,0};
		
		// Use the correct version of getDepth
		if (this.readPosList.size() > 0)
		{
			out = this.readPosList.get(0).getDepth(readPosList, this.apOptions);
		}
		return out;
	}

	
	public  AlignmentPosSnps clone(boolean flip_strand)
	{
		AlignmentPosSnps ap = new AlignmentPosSnps(this.getRefFlipped(), this.chr, this.pos, this.apOptions);

		Vector<ReadPos> newReadPosList = new Vector<ReadPos>(this.readPosList.size());
		Iterator<ReadPos> it = this.readPosList.iterator();
		while (it.hasNext())
		{
			newReadPosList.add(it.next().clone());
		}
		ap.readPosList = newReadPosList;
		
		return ap;
	}
	

	/*****************
	 * 
	 * Setters
	 * 
	 */
	
	
	public void add(Symbol inC, boolean inReadForwardStrand)
	{ 
		ReadPos rp = new ReadPos(inC, inReadForwardStrand);
		this.add(rp);
	}
	
	public void add(ReadPos rp)
	{
		this.readPosList.add(rp.clone());
	}
	
	public void add(Symbol inC, boolean inReadForwardStrand, int inPos, int inQual)
	{
		ReadPos rp = new ReadPosRich(inC, inReadForwardStrand, inPos, inQual);
		this.add(rp);
	}
	
	public void add(ReadPos rp, int inPos, int inQual)
	{
		this.readPosList.add(new ReadPosRich(rp, inPos, inQual));
	}


	
	public void reset()
	{
		this.readPosList = new Vector<ReadPos>();
	}
	


}
