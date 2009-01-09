package edu.usc.epigenome.genomeLibs;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Symbol;

/** 
 * This is a copy of the SymbolCounter class, but in this one we include
 * an arbitrary stratification string to stratify counts
 * @author benb
 *
 */

// Use String name of Symbol as hash key, because Symbol is not iself Comparable
public class SymbolCounterStratified extends TreeMapCounter<String> {

	private static final long serialVersionUID = 8275939009040141042L;

	public SymbolCounterStratified() {
	}

	
	/*** GETTERS ***/
	
	public int getCytosineCount(String strat)
	{
		return getCount(DNATools.c(), strat);
	}
	
	public int getAdenineCount(String strat)
	{
		return getCount(DNATools.a(), strat);
	}
	
	public int getGuanineCount(String strat)
	{
		return getCount(DNATools.g(), strat);
	}
	
	public int getThymineCount(String strat)
	{
		return getCount(DNATools.t(), strat);
	}
	
	public int getCount(Symbol sym, String strat) {
		String key = keyFromSym(sym, strat);
		return super.getCount(key);
	}

	/***
	 * Special function for bisulfite converted DNA.  Returns
	 * @return count(C)/[count(T)+count(C)]
	 */
	public double getRetainedFrac(String strat)
	{
		int c = getCytosineCount(strat);
		int t = getThymineCount(strat);
			
		return ((double)c / ( (double)c + (double)t )); 
	}
	
	/*** SETTERS ***/
	
	public void increment(Symbol sym, String strat, int numToAdd) {
		String key = keyFromSym(sym, strat);
		super.increment(key, numToAdd);
	}

	public void increment(Symbol sym, String strat) {
		String key = keyFromSym(sym, strat);
		super.increment(key);
	}

	public void setCount(Symbol sym, String strat, int newCount) {
		String key = keyFromSym(sym, strat);
		super.setCount(key, newCount);
	}

	/*** STATIC FUNCTION FOR GENERATING HASH KEYS ***/
	protected static String keyFromSym(Symbol sym, String strat)
	{
		String key = sym.getName() + "__" + strat;
		return key;
	}

}
