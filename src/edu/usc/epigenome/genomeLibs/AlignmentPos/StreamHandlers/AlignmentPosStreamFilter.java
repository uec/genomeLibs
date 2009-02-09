/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

/**
 * @author benb
 *
 */
public abstract class AlignmentPosStreamFilter implements
		AlignmentPosStreamHandler {

	/**
	 * @param streamPos
	 * @return true if AP should be passed
	 */

	public abstract boolean elementPasses(AlignmentPosStreamerPosition streamPos);

	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#finish()
	 */
	public void finish() {
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#init()
	 */
	public void init() {
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(edu.usc.epigenome.genomeLibs.AlignmentPos[], edu.usc.epigenome.genomeLibs.AlignmentPos, edu.usc.epigenome.genomeLibs.AlignmentPos[])
	 */
	public boolean streamElement(AlignmentPosStreamerPosition streamPos) {
		return this.elementPasses(streamPos); 
	}

}
