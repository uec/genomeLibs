package edu.usc.epigenome.genomeLibs;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Symbol;

//Use String name of Symbol as hash key, because Symbol is not iself Comparable
public class SymbolCounter extends TreeMapCounter<String> {


	private static final long serialVersionUID = 3083285853454025453L;

	public SymbolCounter() {
	}

	/*** GETTERS ***/
	
	public int getCytosineCount()
	{
		return getCount(DNATools.c());
	}
	
	public int getAdenineCount()
	{
		return getCount(DNATools.a());
	}
	
	public int getGuanineCount()
	{
		return getCount(DNATools.g());
	}
	
	public int getThymineCount()
	{
		return getCount(DNATools.t());
	}
	
	/***
	 * Special function for bisulfite converted DNA.  Returns
	 * @return count(T)/[count(T)+count(C)]
	 */
	public double getConvertedFrac()
	{
		int c = getCytosineCount();
		int t = getThymineCount();
			
		return ((double)t / ( (double)c + (double)t )); 
	}
	
	public int getCount(Symbol sym) {
		String key = keyFromSym(sym);
		return super.getCount(key);
	}

	/*** SETTERS ***/
	
	public void increment(Symbol sym, int numToAdd) {
		String key = keyFromSym(sym);
		super.increment(key, numToAdd);
	}

	public void increment(Symbol sym) {
		String key = keyFromSym(sym);
		super.increment(key);
	}

	public void setCount(Symbol sym, int newCount) {
		String key = keyFromSym(sym);
		super.setCount(key, newCount);
	}
	
	protected static String keyFromSym(Symbol sym)
	{
		String key = sym.getName();
		return key;
	}
	

}
