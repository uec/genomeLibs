/**
 * 
 */
package edu.usc.epigenome.genomeLibs;

import java.util.*;

/**
 * @author benb
 *
 */
public interface AlignmentPosStreamHandler {

	/** Constructors **/
	
	
	/**
	 * Called before any elements are streamed 
	 */
	public abstract void init();
	
		
	/**
	 * Called after all elements are streamed
	 */
	public abstract void finish();
	
	/**
	 * @param priorAps is a list APs preceeding current one (length determined by streamer)
	 * @param currentAp is the current AP
	 * @param nextAps is a list APs following current one (length determined by streamer)
	 * @return true if AP should be passed to next handler
	 */
	public abstract boolean streamElement(LinkedList<AlignmentPos> priorAps, 
			AlignmentPos currentAp, LinkedList<AlignmentPos> nextAps);
	
	
}
