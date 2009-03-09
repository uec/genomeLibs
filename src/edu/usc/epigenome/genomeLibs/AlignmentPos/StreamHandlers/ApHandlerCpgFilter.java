package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import edu.usc.epigenome.genomeLibs.CpgPair;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

abstract public class ApHandlerCpgFilter extends APHandlerCpgHandler {

	public ApHandlerCpgFilter() {
	}

	@Override
	public boolean streamCpgPair(AlignmentPosStreamerPosition streamPos,
			CpgPair pair) {
		fReportCounts = false;
		return this.elementPasses(streamPos, pair);
	}
	
	
	
	abstract public boolean elementPasses(AlignmentPosStreamerPosition streamPos,CpgPair pair);

	
}
