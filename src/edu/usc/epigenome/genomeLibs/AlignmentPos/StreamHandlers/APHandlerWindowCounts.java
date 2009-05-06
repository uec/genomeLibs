/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.GregorianCalendar;

import edu.usc.epigenome.genomeLibs.AlignmentPos.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;
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
	public APHandlerWindowCounts(int inWindSize, boolean strandSpecific, String inGenome) {
		super(inWindSize, strandSpecific, inGenome);
	}

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	public void init() {
	}

	public void finish() {
		System.err.println((new GregorianCalendar()).getTime() + "\tAPHandlerWindowCounts finished");
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	public boolean streamElement(AlignmentPosStreamerPosition streamPos) 
	{
		this.increment(streamPos.currentAp);
		return true;
	}


}
