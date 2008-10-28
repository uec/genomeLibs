/**
 * 
 */
package edu.usc.epigenome.genomeLibs;

import java.util.*;

/**
 * @author benb
 *
 */
public abstract class AlignmentPosStreamFilter extends
		AlignmentPosStreamHandler {

	/**
	 * 
	 */
	public AlignmentPosStreamFilter() {
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#finish()
	 */
	@Override
	public void finish() {

	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#init()
	 */
	@Override
	public void init() {
	}

	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	@Override
	public boolean streamElement(LinkedList<AlignmentPos> priorAps,
			AlignmentPos currentAp, LinkedList<AlignmentPos> nextAps) {
		return this.elementPasses(priorAps, currentAp, nextAps);
	}

	/**
	 * @param priorAps is a list APs preceeding current one (length determined by streamer)
	 * @param currentAp is the current AP
	 * @param nextAps is a list APs following current one (length determined by streamer)
	 * @return true if AP should be passed
	 */

	abstract protected boolean elementPasses(LinkedList<AlignmentPos> priorAps, 
			AlignmentPos currentAp, LinkedList<AlignmentPos> nextAps);

}
