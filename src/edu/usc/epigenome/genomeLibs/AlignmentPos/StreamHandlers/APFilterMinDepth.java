package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import org.biojava.bio.seq.DNATools;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

public class APFilterMinDepth extends AlignmentPosStreamFilter {

	protected int minDepth = 1;
	protected boolean eachStrand = false;
	
	public APFilterMinDepth(int inMinDepth, boolean inEachStrand) {
		minDepth = inMinDepth;
		eachStrand = inEachStrand;
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



	/**
	 * @return the eachStrand
	 */
	public boolean isEachStrand() {
		return eachStrand;
	}

	/**
	 * @param eachStrand the eachStrand to set
	 */
	public void setEachStrand(boolean eachStrand) {
		this.eachStrand = eachStrand;
	}

	@Override
	public boolean elementPasses(AlignmentPosStreamerPosition streamPos) {

		boolean passes = true;
		if (this.isEachStrand())
		{
			passes &= (streamPos.currentAp.getDepth(true) >= this.getMinDepth()); 
			passes &= (streamPos.currentAp.getDepth(false) >= this.getMinDepth()); 
		}
		else
		{
			passes = (streamPos.currentAp.getTotalDepth() >= this.getMinDepth());
		}
		
		return passes;
	}
}
