package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import org.biojava.bio.seq.DNATools;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

public class APFilterNmer extends AlignmentPosStreamFilter {

	protected String nmer = null;
	protected boolean fwStrandOnly = true;
	
	public APFilterNmer(String inNmer, boolean inFwStrandOnly) {
		nmer = inNmer;
		fwStrandOnly = inFwStrandOnly;
	}

	@Override
	public boolean elementPasses(AlignmentPosStreamerPosition streamPos) {

		int n = nmer.length();
		
		boolean passes = false;
		if (!fwStrandOnly)
		{
			System.err.println("APFilterNmer does not yet implement fwStrandOnly==0");
			System.exit(0);
		}
		else
		{
			if (streamPos.nextAps.length < (n-1))
			{
				System.err.println("APFilterNmer must have at least 2 bp of sequence context on either side");
				//	(new Exception()).printStackTrace();
				System.exit(0);
			}
			
			// Is it the nmer
			String actual = streamPos.currentAp.getRefToken() + AlignmentPos.getRefTokens(streamPos.nextAps).substring(0, n-1);
			//System.err.println(actual + " == " + nmer + " ?");
			
			passes = (actual.equalsIgnoreCase(nmer));
			
		}

		return passes;
	}
}
