/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.Counters.TreeMapCounter;


/**
 * @author benb
 * 
 *
 */
public class APHandlerDepthCounts extends TreeMapCounter<Integer> implements AlignmentPosStreamHandler {


	private static final long serialVersionUID = 9000604824830397400L;

	/**
	 * Constructor
	 */
	public APHandlerDepthCounts() {
		super();
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
		// Check how many reads overlap this position, do forward and reverse separately.
		int fw = currentAp.getDepth(true);
		int rev = currentAp.getDepth(false);
		
		if (fw>0)
		{
			this.increment(new Integer(fw));
		}
		
		if (rev>0)
		{
			// Make reverse negative
			this.increment(new Integer(-rev));
		}
		return true;
	}

	
	

}
