package edu.usc.epigenome.genomeLibs.SampleOrdering;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;

import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleAlphabet;
import org.biojava.bio.symbol.Symbol;

import edu.usc.epigenome.genomeLibs.MiscUtils;

public class SampleOrderingAlphabet extends SimpleAlphabet {

	// So we can get numerics for plotting on a wig track
	protected Map<String,Integer> symbolNums = new HashMap<String,Integer>(500);
	
	public SampleOrderingAlphabet(String name, Set<Character> sampleIds)
	throws Exception
	{
		super(name);
		
		Set<String> symStrs = SampleOrderingAlphabet.generateAllSymbolStrings(sampleIds);
		
		int onNum = 1;
		for (String symStr : symStrs)
		{
			Symbol sym = AlphabetManager.createSymbol(symStr);
			this.addSymbol(sym);
			
			symbolNums.put(symStr, onNum++);
		}
	}

	
	public Symbol getSymbol(char[] sampleIds, double[] vals, double slopForTies)
	{
		String str = getSymbolString(sampleIds, vals, slopForTies);
		Symbol sym = null;
		try
		{
			sym = this.getTokenization("name").parseToken(str);
		}
		catch (Exception e)
		{
			System.err.println("Why can't we find symbol \"" + sym + "\"?");
			e.printStackTrace();
			System.exit(1);
		}
		return sym;
	}

	public int getSymbolNumber(char[] sampleIds, double[] vals, double slopForTies)
	{
		String str = getSymbolString(sampleIds, vals, slopForTies);
		//System.err.println("Got string : " + str);
		return symbolNums.get(str);
	}
	
	public static String getSymbolString(char[] sampleIds, double[] vals, double slopForTies)
	{
		// Remember that in ties, the smaller letter has to go first, i.e. "abc>d" not "acb>d" 
		
		int nV = vals.length; 
		
		// Is there a faster way to get sort order?
		List<IdValPair> pairs = new ArrayList<IdValPair>(nV);
		for (int i = 0; i < nV; i++)
		{
			IdValPair pair = new IdValPair(sampleIds[i],vals[i]);
			pairs.add(pair);
		}
		
		Collections.sort(pairs, new SampleOrderingAlphabet.ValComparator());
		
		StringBuffer outBuf = new StringBuffer(nV*2);
		
		// We have to keep track of pairs so we can sort them by id.
		List<IdValPair> tiedPairs = new ArrayList<IdValPair>(nV);
		tiedPairs.add(pairs.get(0));
		for (int i = 1; i < nV; i++)
		{
			IdValPair lastPair = pairs.get(i-1);
			IdValPair pair = pairs.get(i);
			
			boolean tie = (pair.val + slopForTies) >= lastPair.val;
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(
					String.format("Comparing %s to %s, tie=%s out=%s\n", lastPair, pair, ""+tie, outBuf.toString()));
			
			if (!tie)
			{
				// Output the previous tie group and clear it for the
				// next go round.
				Collections.sort(tiedPairs, new SampleOrderingAlphabet.IdComparator());
				for (IdValPair outPair : tiedPairs)
				{
					outBuf.append(outPair.id);
				}
				tiedPairs.clear();

				// Output the GT unless this is the last one
				if (i<nV) outBuf.append(">");
			}

			// Add ourself.
			tiedPairs.add(pair);
		}
		
		// Add the last tied group
		Collections.sort(tiedPairs, new SampleOrderingAlphabet.IdComparator());
		for (IdValPair outPair : tiedPairs)
		{
			outBuf.append(outPair.id);
		}

		
		return(outBuf.toString());
	}
	
	protected static Set<String> generateAllSymbolStrings(Set<Character>sampleIds)
	{
		Set<String> out = new HashSet<String>(3*(int)MiscUtils.factorial(sampleIds.size()));
		generateSymbolStrings(out, "", sampleIds);
		return out;
	}
	
	// Adds final leaf nodes to out.  Better than passing big lists around
	protected static void generateSymbolStrings(Set<String> out, String prefix, Set<Character> remainingIds)
	{
		int nrem = remainingIds.size();
		
		if (nrem==0)
		{
			// Done
			out.add(prefix);
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).finer("Adding to alphabet: " + prefix);
			return;
		}
		
		int preLen = prefix.length();
		boolean priorWasGt = ((preLen>0) && (prefix.charAt(preLen-1) == '>'));
		
		if ((preLen>0) && !priorWasGt)
		{
			Set<Character> cloneIds = new HashSet(remainingIds);
			generateSymbolStrings(out, prefix + ">", cloneIds);
		}
		
		for (Character c : remainingIds)
		{
			Set<Character> cloneIds = new HashSet(remainingIds);
			cloneIds.remove(c);
			
			// If it's a tie, ordering doesn't matter.  So we will only add ones that 
			// are naturally smaller than the one before.
			if (priorWasGt || (preLen==0) || (prefix.charAt(preLen-1) < c.charValue()))
			{
				generateSymbolStrings(out, prefix + c, cloneIds);
			}
		}
	}
	
	
	public static class IdValPair
	{
		public char id;
		public double val;

		public IdValPair(char id, double val) {
			super();
			this.id = id;
			this.val = val;
		}

		@Override
		public String toString() {
			return String.format("%c.%f", id,val);
		}
	}
	
	public static class IdComparator implements Comparator<IdValPair>
	{
		@Override
		public int compare(IdValPair o1, IdValPair o2) {
			//System.err.println("Sorting ids: " + o1 + ", " + o2);
			if (o1.id > o2.id) return 1; 
			if (o1.id < o2.id) return -1; 
			return 0;
		}
	}
	
	// DESCENDING
	public static class ValComparator implements Comparator<IdValPair>
	{
		@Override
		public int compare(IdValPair o1, IdValPair o2) {
			if (o1.val > o2.val) return -1; 
			if (o1.val < o2.val) return 1; 
			return 0;
		}
	}

}
