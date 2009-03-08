package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import edu.usc.epigenome.genomeLibs.CpgPair;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

public class APHandlerCpgFilterNonSnpCpgs extends ApHandlerCpgFilter {

	public APHandlerCpgFilterNonSnpCpgs() {
	}

	@Override
	public boolean elementPasses(AlignmentPosStreamerPosition streamPos,
			CpgPair pair) {

		double fw = pair.getSnpProb(true);
		double rev = pair.getSnpProb(false);
		
		return ((fw<0.1) && (rev<0.1));
	}

}
