package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import org.biojava.bio.seq.DNATools;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

public class APFilterCytosines extends AlignmentPosStreamFilter {

	public APFilterCytosines() {
		// TODO Auto-generated constructor stub
	}

	@Override
	public boolean elementPasses(AlignmentPosStreamerPosition streamPos) {

		return streamPos.currentAp.getRef().equals(DNATools.c());
	}
}
