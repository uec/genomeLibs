package edu.usc.epigenome.genomeLibs.Counters;

import java.util.TreeMap;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Symbol;

//Use String name of Symbol as hash key, because Symbol is not iself Comparable
public class NmerCounter extends TreeMapCounter<String> {


	private static final long serialVersionUID = -1497021469332805936L;

	// Make a static symbol buffer so we don't always have to reallocate
	private static final int MAX_NMER = 10;
	private static TreeMap<Integer,Symbol[]> sSymSeqBuffers = null;
	
	public NmerCounter() {
		// Init the seq buffers
		sSymSeqBuffers = new TreeMap<Integer,Symbol[]>();
		for (int i = 1; i <= MAX_NMER; i++)
		{
			sSymSeqBuffers.put(new Integer(i), new Symbol[i]);
		}
	}

	/*** GETTERS ***/
	
	public int getCount(Symbol[] syms)
	{
		String key = keyFromSyms(syms);
		return super.getCount(key);
	}
	
	/*** SETTERS ***/
	
	public void incrementAllSubstrings(Symbol[] syms, int maxNmer)
	{
		for (int i = 0; i < syms.length; i++)
		{
			int remaining = syms.length - i;
			
			for (int j = 1; j <= maxNmer; j++)
			{
				if (remaining >= j)
				{
					Symbol[] subSyms = sSymSeqBuffers.get(new Integer(j));
					for (int k = 1; k<=j; k++)
					{
						subSyms[k-1] = syms[i+k-1];
					}
					this.increment(subSyms);
				}
			}
		}
	}
	
	public void increment(Symbol[] syms, int numToAdd) {
		String key = keyFromSyms(syms);
		super.increment(key, numToAdd);
	}

	public void decrement(Symbol[] syms, int numToSubtract) {
		String key = keyFromSyms(syms);
		super.decrement(key, numToSubtract);
	}

	public void increment(Symbol[] syms) {
		String key = keyFromSyms(syms);
		super.increment(key);
	}

	public void decrement(Symbol[] syms) {
		String key = keyFromSyms(syms);
		super.decrement(key);
	}

	public void setCount(Symbol[] syms, int newCount) {
		String key = keyFromSyms(syms);
		super.setCount(key, newCount);
	}
	
	protected static String keyFromSyms(Symbol[] syms)
	{
		StringBuffer buf = new StringBuffer(syms.length);
		try
		{
			for (int i = 0; i < syms.length; i++)
			{
				buf.append(DNATools.dnaToken(syms[i]));
			}
		}
		catch (Exception e)
		{
			System.err.println("NmerCounter called with non-DNA symbol");
			System.exit(1);
		}
		
		return buf.toString();
	}
	

}
