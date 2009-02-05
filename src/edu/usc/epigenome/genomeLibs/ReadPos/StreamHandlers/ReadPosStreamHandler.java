/**
 * 
 */
package edu.usc.epigenome.genomeLibs.ReadPos.StreamHandlers;

import java.util.*;

import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;

/**
 * @author benb
 *
 */
public interface ReadPosStreamHandler {

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
	 * @param currentRp is the current RP
	 * @return true if AP should be passed to next handler
	 */
	public abstract boolean streamElement(ReadPos currentAp);
	
	
}
