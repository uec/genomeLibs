package edu.usc.epigenome.genomeLibs;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

public class APFilterCphs extends AlignmentPosStreamFilter {

	public static SimpleAlphabet alph;
	
	// Initializae alph
	static
	{
		try
		{
		alph = new SimpleAlphabet();
		alph.addSymbol(DNATools.a());
		alph.addSymbol(DNATools.c());
		alph.addSymbol(DNATools.t());
		}
		catch (IllegalSymbolException e)
		{
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	public APFilterCphs() {
	}

	@Override
	public boolean elementPasses(AlignmentPos[] priorAps,
			AlignmentPos currentAp, AlignmentPos[] nextAps) {

		if ((nextAps.length < 1) || (priorAps.length < 1))
		{
			System.err.println("APFilterCphs must have at least bp of sequence context on either side");
		//	(new Exception()).printStackTrace();
			System.exit(0);
		}
		
		
		return ((currentAp.getRef().equals(DNATools.c())) &&
				(alph.contains(nextAps[0].getRef())));
	}
}
                                                        