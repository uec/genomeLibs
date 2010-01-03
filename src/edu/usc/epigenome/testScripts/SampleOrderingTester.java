package edu.usc.epigenome.testScripts;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.biojava.bio.symbol.Symbol;


import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.SampleOrdering.SampleOrderingAlphabet;

public class SampleOrderingTester {

	/**
	 * @param args
	 */
	public static void main(String[] args) 
	throws Exception {
		// TODO Auto-generated method stub

		Set<Character> set = new HashSet(Arrays.asList('a','b','c','d'));
		SampleOrderingAlphabet alph = new SampleOrderingAlphabet("ab list", set);
		Iterator syms = alph.iterator();
		System.err.print("Alph ("+ alph.size() + ") =");
		while (syms.hasNext())
		{
			Symbol sym = (Symbol)syms.next();
			System.err.print(", " + sym.getName());
		}
		System.err.println("");
		char[] ids = {'a','b','c','d'};
		double[] vals = {10.0, 9.9, 16.0, 40.0};
		System.err.println("Symbol for abc>d:\n " + alph.getSymbol(ids, vals, 0.0));
	}

}
