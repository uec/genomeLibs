/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

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
	 * 
	 * priorAps and nextAps are guaranteed to contain APs which are adjacent to each other 
	 * on the chromosome.
	 */
	public abstract boolean streamElement(AlignmentPosStreamerPosition streamPos);
	
	
}
