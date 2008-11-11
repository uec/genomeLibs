/**
 * 
 */
package edu.usc.epigenome.genomeLibs;

import java.util.*;

/**
 * @author benb
 *
 */
public interface AlignmentPosStreamFilter extends
		AlignmentPosStreamHandler {

	/**
	 * @param priorAps is a list APs preceeding current one (length determined by streamer)
	 * @param currentAp is the current AP
	 * @param nextAps is a list APs following current one (length determined by streamer)
	 * @return true if AP should be passed
	 */

	public abstract boolean elementPasses(LinkedList<AlignmentPos> priorAps, 
			AlignmentPos currentAp, LinkedList<AlignmentPos> nextAps);

}
