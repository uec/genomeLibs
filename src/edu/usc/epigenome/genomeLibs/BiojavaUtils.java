package edu.usc.epigenome.genomeLibs;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.StrandedFeature;
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

	public static int strandToInt(StrandedFeature.Strand strand)
	{
		int out = 0;
		
		if (strand == StrandedFeature.NEGATIVE)
		{
			out = -1;
		}
		else if (strand == StrandedFeature.POSITIVE)
		{
			out = 1;
		}
		
		return out;
	}
	
	
}
