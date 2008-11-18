/**
 * 
 */
package edu.usc.epigenome.genomeLibs;

/**
 * @author benb
 *
 */
public abstract class AlignmentPosStreamFilter implements
		AlignmentPosStreamHandler {

	/**
	 * @param priorAps is a list APs preceeding current one (length determined by streamer)
	 * @param currentAp is the current AP
	 * @param nextAps is a list APs following current one (length determined by streamer)
	 * @return true if AP should be passed
	 */

	public abstract boolean elementPasses(AlignmentPos[] priorAps, 
			AlignmentPos currentAp, AlignmentPos[] nextAps);

	
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
	public boolean streamElement(AlignmentPos[] priorAps,
			AlignmentPos currentAp, AlignmentPos[] nextAps) {
		return this.elementPasses(priorAps, currentAp, nextAps); 
	}

}
