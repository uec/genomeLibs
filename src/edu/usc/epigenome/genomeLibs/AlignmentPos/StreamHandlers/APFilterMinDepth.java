package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import org.biojava.bio.seq.DNATools;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;

public class APFilterMinDepth extends AlignmentPosStreamFilter {

	protected int minDepth = 1;
	
	public APFilterMinDepth(int inMinDepth) {
		minDepth = inMinDepth;
	}

	
	
	/**
	 * @return the minDepth
	 */
	public int getMinDepth() {
		return minDepth;
	}



	/**
	 * @param minDepth the minDepth to set
	 */
	public void setMinDepth(int minDepth) {
		this.minDepth = minDepth;
	}



	@Override
	public boolean elementPasses(AlignmentPos[] priorAps,
			AlignmentPos currentAp, AlignmentPos[] nextAps) {

		return (currentAp.getTotalDepth() >= this.getMinDepth());
	}
}
