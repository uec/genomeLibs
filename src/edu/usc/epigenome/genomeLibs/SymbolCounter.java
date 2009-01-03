package edu.usc.epigenome.genomeLibs;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Symbol;

// Use Character because underlying Symbol is not Comparable
public class SymbolCounter extends TreeMapCounter<String> {

	private static final long serialVersionUID = 8275939009040141042L;

	public SymbolCounter() {
	}

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
	 * #T/(#C+#T)
	 * @return
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
	
//	protected static Character keyFromSym(Symbol sym)
//	{
//		Character key = new Character(BiojavaUtils.dnaTokenNoException(sym));
//		return key;
//	}

}
