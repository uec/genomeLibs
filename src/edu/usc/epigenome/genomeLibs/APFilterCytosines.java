package edu.usc.epigenome.genomeLibs;

import org.biojava.bio.seq.DNATools;

public class APFilterCytosines extends AlignmentPosStreamFilter {

	public APFilterCytosines() {
		// TODO Auto-generated constructor stub
	}

	@Override
	public boolean elementPasses(AlignmentPos[] priorAps,
			AlignmentPos currentAp, AlignmentPos[] nextAps) {

		return currentAp.getRef().equals(DNATools.c());
	}
}
