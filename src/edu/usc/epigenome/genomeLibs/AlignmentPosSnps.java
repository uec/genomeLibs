package edu.usc.epigenome.genomeLibs;

import java.util.*;

import org.biojava.bio.seq.StrandedFeature;
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

	
	
	/**
	 * @return a list of read positions , with identical ones removed
	 * according to apOptions.maxIdentical
	 */
	public Vector<ReadPos> getReadPositions()
	{
		
		TreeMap<String,Integer> counts = new TreeMap<String,Integer>();
		Iterator<ReadPos> it = readPosList.iterator();
		Vector<ReadPos> out = new Vector<ReadPos>(readPosList.size());
		while (it.hasNext())
		{
			ReadPos rp = it.next();
			boolean add = true;
			
			if (this.apOptions.maxIdentical > 0)
			{
				
				int cycle = rp.getCycle();
				if (cycle != ReadPos.UNKNOWN_CYCLE) // If it's unknown , can't eliminate
				{
					String key = rp.getStrand() + "__" + cycle;
					//System.err.println("Looking for key: " + key);
					int val = (counts.get(key) == null) ? 0 : ((Integer)counts.get(key)).intValue();
					val++; // The current one
					add = (val <= this.apOptions.maxIdentical);
					counts.put(key, new Integer(val));
				}
			}
			
			
			if (add)
			{
				out.add(rp);
			}
		}
			
		return out;
	}
	
	
	
	
	
	public int[] getDepth()
	{
		int [] out = new int[] {0,0};
		
		Iterator<ReadPos> rpIt = this.getReadPositions().iterator();
		while (rpIt.hasNext())
		{
			ReadPos rp = rpIt.next();
			int index = (rp.getStrand() == StrandedFeature.NEGATIVE) ? 1 : 0;
			out[index]++;
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
