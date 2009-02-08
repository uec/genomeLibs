/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import edu.usc.epigenome.genomeLibs.AlignmentPos.*;

/**
 * @author benb
 * 
 * You pass in a filterer and a hanlder, and this object will only
 * process APs passing the filterer.
 *
 */
public class AlignmentPosStreamFilteredHandler implements
		AlignmentPosStreamHandler {

	protected AlignmentPosStreamFilter filterer = null;
	protected AlignmentPosStreamHandler handler = null;
	
	
	
	/**
	 * @param filterer
	 * @param handler
	 */
	public AlignmentPosStreamFilteredHandler(AlignmentPosStreamFilter filterer,
			AlignmentPosStreamHandler handler) {
		super();
		this.filterer = filterer;
		this.handler = handler;
	}

	
	
	
	

	/**
	 * @return the handler
	 */
	public AlignmentPosStreamHandler getHandler() {
		return handler;
	}






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
		boolean out = filterer.elementPasses(priorAps, currentAp, nextAps);
		if (out) handler.streamElement(priorAps, currentAp, nextAps);
		return true;
	}

}
