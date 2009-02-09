/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import edu.usc.epigenome.genomeLibs.AlignmentPos.*;
import edu.usc.epigenome.genomeLibs.Counters.*;


/**
 * @author benb
 * 
 *
 */
public class APHandlerWindowCounts extends AlignmentPosWindCounter implements AlignmentPosStreamHandler {

	private static final long serialVersionUID = 8989385550284144139L;

	/**
	 * Constructor
	 */
	public APHandlerWindowCounts(int inWindSize, boolean strandSpecific) {
		super(inWindSize, strandSpecific);
	}

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	public void init() {
	}

	public void finish() {
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	public boolean streamElement(AlignmentPos[] pre, AlignmentPos currentAp, AlignmentPos[] post) 
	{
		this.increment(currentAp);
		return true;
	}


}
