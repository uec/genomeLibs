/**
 * 
 */
package edu.usc.epigenome.genomeLibs;

import java.util.Collection;
import java.util.Iterator;
import java.util.TreeMap;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Symbol;

/**
 * @author benb
 *
 */
public class ReadPosRich extends ReadPos {

	/* Obj vars */
	protected int pos = UNKNOWN;
	protected int qual = UNKNOWN;
	
	/** Constructors **/
	
	public ReadPosRich(Symbol inSym, boolean inForwardStrand, int inPos, int inQual)
	{
		this.sym = inSym;
		this.strand = (inForwardStrand) ? StrandedFeature.POSITIVE : StrandedFeature.NEGATIVE;
		this.qual = inQual;
		this.pos = inPos;
	}

	public ReadPosRich(ReadPos inRp, int inPos, int inQual)
	{
		this.sym = inRp.getSym();
		this.strand = inRp.getStrand();
		this.qual = inQual;
		this.pos = inPos;
	}
	
	
	/* Getters/Setters */

	/**
	 * @return the readPos
	 */
	public int getReadPos() {
		return pos;
	}


	/**
	 * @param readPos the readPos to set
	 */
	public void setReadPos(int pos) {
		this.pos = pos;
	}


	/**
	 * @return the qual
	 */
	public int getQual() {
		return qual;
	}


	/**
	 * @param qual the qual to set
	 */
	public void setQual(int qual) {
		this.qual = qual;
	}
	
	
	
	
	
	
	/**
	 * @param a list of ReadPos objs
	 * @return int[0]=fw_depth , int[1]=rev_depth
	 */

	public int[] getDepth(Collection<ReadPos> posList, ReadPosOptions ro)
	{
		int[] depth = new int[] {0,0};

		TreeMap<String,Integer> counts = new TreeMap<String,Integer>();
		
		Iterator<ReadPos> it = posList.iterator();
		while (it.hasNext())
		{
			ReadPosRich rp = (ReadPosRich)it.next();
			boolean add = true;
			
			if (ro.maxIdentical > 0)
			{
				String key = rp.getStrand() + "__" + rp.getReadPos();
				//System.err.println("Looking for key: " + key);
				int val = (counts.get(key) == null) ? 0 : ((Integer)counts.get(key)).intValue();
				val++; // The current one
				add = (val <= ro.maxIdentical);
				counts.put(key, new Integer(val));
			}
			
			
			if (add)
			{
				int index = (rp.getStrand()==StrandedFeature.NEGATIVE) ? 1 : 0;
				depth[index]++;
			}
		}
		
		return depth;
	}
	

	
	
	
	
}
