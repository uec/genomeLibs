package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import org.biojava.bio.seq.DNATools;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

public class APFilterFixedIntervals extends AlignmentPosStreamFilter {

	protected int interval = 1;
	protected boolean eachStrand = false;
	
	public APFilterFixedIntervals(int inInterval) {
		interval = inInterval;
	}

	/**
	 * @return the minDepth
	 */
	public int getInterval() {
		return interval;
	}



	/**
	 * @param minDepth the minDepth to set
	 */
	public void setInterval(int inInterval) {
		this.interval = inInterval;
	}



	

	@Override
	public boolean elementPasses(AlignmentPosStreamerPosition streamPos) {

		boolean passes = true;

		int pos = streamPos.currentAp.getPos();
		passes = ((pos % interval) == 0);
		//System.err.println(pos + " % " + interval + " , passes= " + passes);
		
		return passes;
	}
}
