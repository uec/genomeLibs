package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

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
	public boolean elementPasses(AlignmentPosStreamerPosition streamPos) {

		if ((streamPos.nextAps.length < 1) || (streamPos.priorAps.length < 1))
		{
			System.err.println("APFilterCphs must have at least bp of sequence context on either side");
		//	(new Exception()).printStackTrace();
			System.exit(0);
		}
		
		
		return ((streamPos.currentAp.getRef().equals(DNATools.c())) &&
				(alph.contains(streamPos.nextAps[0].getRef())));
	}
}
                                                        