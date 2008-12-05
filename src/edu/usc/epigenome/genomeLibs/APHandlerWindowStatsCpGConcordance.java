package edu.usc.epigenome.genomeLibs;

import java.util.*;

import org.biojava.bio.seq.StrandedFeature;

public class APHandlerWindowStatsCpGConcordance extends APHandlerWindowStats {

	public APHandlerWindowStatsCpGConcordance(int inWindSize) {
		super(inWindSize);
	}

	@Override
	public void init() {
		super.init();
	}

	@Override
	public void finish() {
		super.finish();
	}


	@Override
	public boolean streamWindow(AlignmentPos[] priorAps,
			AlignmentPos currentAp, AlignmentPos[] nextAps,
			Queue<AlignmentPos> apWind) {

		boolean currentFw = (currentAp.getStrand() == StrandedFeature.POSITIVE);
		LinkedList<AlignmentPos> apWindCopy = new LinkedList<AlignmentPos>(); // Clone
		apWindCopy.addAll(apWind);
		AlignmentPos.removeApsByStrand(apWindCopy, !currentFw);
		
		System.err.println(
				currentAp + 
				"\t" + 
				AlignmentPos.getRefTokens(priorAps) + 
				"/" + 
				currentAp.getRefToken() + 
				"/" +
				AlignmentPos.getRefTokens(nextAps) + 
				"\t" + 
				AlignmentPos.getRefTokens(apWindCopy) +
				"\t" + 
				AlignmentPos.getPosString(apWindCopy)
				);
		
		return true;
	}

}
