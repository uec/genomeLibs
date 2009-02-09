package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import org.biojava.bio.seq.DNATools;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

public class APFilterCpgs extends AlignmentPosStreamFilter {

	public APFilterCpgs() {
		// TODO Auto-generated constructor stub
	}

	@Override
	public boolean elementPasses(AlignmentPosStreamerPosition streamPos) {

		if ((streamPos.nextAps.length < 1) || (streamPos.priorAps.length < 1))
		{
			System.err.println("APFilterCpgs must have at least bp of sequence context on either side");
		//	(new Exception()).printStackTrace();
			System.exit(0);
		}
		
		return ((streamPos.currentAp.getRef().equals(DNATools.c())) &&
				(streamPos.nextAps[0].getRef().equals(DNATools.g())));
	}
}
