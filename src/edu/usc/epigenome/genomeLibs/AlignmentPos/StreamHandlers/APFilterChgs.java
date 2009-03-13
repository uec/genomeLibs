package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import org.biojava.bio.seq.DNATools;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

public class APFilterChgs extends AlignmentPosStreamFilter {

	public APFilterChgs() {
		// TODO Auto-generated constructor stub
	}

	@Override
	public boolean elementPasses(AlignmentPosStreamerPosition streamPos) {

		if ((streamPos.nextAps.length < 2) || (streamPos.priorAps.length < 2))
		{
			System.err.println("APFilterChgs must have at least 2 bp of sequence context on either side");
		//	(new Exception()).printStackTrace();
			System.exit(0);
		}
		
		return ((streamPos.currentAp.getRef().equals(DNATools.c())) &&
				(!streamPos.nextAps[0].getRef().equals(DNATools.g())) &&
				(streamPos.nextAps[1].getRef().equals(DNATools.g())));	
	}
}
