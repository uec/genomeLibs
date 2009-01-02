package edu.usc.epigenome.genomeLibs;

import java.util.*;

import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.*;

public class AlignmentPosSnps extends AlignmentPos {

	/* Obj vars */
	protected Vector<ReadPos> readPosList = new Vector<ReadPos>();
	


	/*****************
	 *  Constructors
	 */

	public AlignmentPosSnps() {
		super();
	}

	
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

	public Vector<ReadPos> getReadPositions(boolean fwOnly)
	{
		
		TreeMap<String,Integer> counts = new TreeMap<String,Integer>();
		Iterator<ReadPos> it = readPosList.iterator();
		Vector<ReadPos> out = new Vector<ReadPos>(readPosList.size());
		while (it.hasNext())
		{
			ReadPos rp = it.next();
			boolean add = true;

			if (fwOnly && (rp.strand == StrandedFeature.NEGATIVE))
			{

			}
			else
			{

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

	
	public AlignmentPos clone(boolean flipStrand)
	{
		// AlignmentPosSnps ap = new AlignmentPosSnps(this.getRef(!flipStrand), this.chr, this.pos, this.apOptions);
		
		// Create a new instance of the correct type
		AlignmentPosSnps ap;
		try {
			Class<? extends AlignmentPosSnps> c = this.getClass();
			ap = c.newInstance();
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		} 
		ap.setRef(this.getRef(!flipStrand));
		ap.setChr(this.getChr());
		ap.setPos(this.getPos());
		ap.setApOptions(this.getApOptions());
		
		// Add the RPs
		Vector<ReadPos> newReadPosList = new Vector<ReadPos>(this.readPosList.size());
		Iterator<ReadPos> it = this.readPosList.iterator();
		while (it.hasNext())
		{
			ReadPos rp = it.next();
			if (flipStrand) rp = rp.reverseComplement();
			newReadPosList.add(rp);
		}
		ap.readPosList = newReadPosList;
		
		// And change the strand if necessary
		StrandedFeature.Strand strand = (flipStrand) ? this.getStrand().flip() : this.getStrand();
		ap.setStrand(strand);
//		System.err.println("Cloning. old strand=" + this.getStrand() + "\tnew=" + strand);
//		System.err.println(this.toString());
//		System.err.println(ap.toString());
		
		return ap;
	}
	


	/*****************
	 * 
	 * Setters
	 * 
	 */
	
	
	public void add(Symbol inC, StrandedFeature.Strand inStrand)
	{ 
		ReadPos rp = new ReadPos(inC, inStrand);
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


	public void removeRevStrandReads()
	{
		Vector<ReadPos> newRp = this.getReadPositions(true);
		this.readPosList = newRp;
	}
	
	public void reset()
	{
		this.readPosList = new Vector<ReadPos>();
	}



}
