/**
 * 
 */
package edu.usc.epigenome.genomeLibs;

import java.util.*;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

/**
 * @author benb
 * 
 *
 */
public class RPHandlerSymbolCounts extends ReadPosCounts implements ReadPosStreamHandler {

	/**
	 * 
	 */
	public RPHandlerSymbolCounts() {
	}

	

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	public boolean streamElement(ReadPos currentRp) 
	{
		boolean passes = true;

		//System.err.println("Streaming element: " + currentRp.commaSeparatedLine());
			Integer count = cycleCounts.get(currentRp);
			if (count == null) count = new Integer(0);
			count++;
			cycleCounts.put(currentRp, count);
		
		return passes;
	}

	
	
	

}
