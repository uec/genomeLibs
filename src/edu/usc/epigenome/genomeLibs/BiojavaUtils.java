package edu.usc.epigenome.genomeLibs;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Symbol;

public class BiojavaUtils {


	public static char dnaTokenNoException(Symbol sym)
	{
		char c = '0';
		try
		{
			c = DNATools.dnaToken(sym);
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		
		return c;	
	}


	
}
