/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

/**
 * @author benb
 * 
 *
 */
abstract public class APHandlerWindowStats implements AlignmentPosStreamHandler {

	public LinkedList<AlignmentPos> window;
	public int windSize = 0;
	
	/**
	 * 
	 */
	public APHandlerWindowStats(int inWindSize) 
	{
		windSize = inWindSize;
	}


	public void init()
	{
		window = new LinkedList<AlignmentPos>();
	}


	public void finish()
	{
		window = null;
	}



	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	public boolean streamElement(AlignmentPosStreamerPosition streamPos) 
	{

		// Eat APs off the head until we are within range.
		String curChr = streamPos.currentAp.getChr();
		int curPos = streamPos.currentAp.getPos();
		boolean done = false;
		AlignmentPos endAp;
		while (!done && ((endAp = window.peek()) != null))
		{
			if ((endAp.getChr().equals(curChr)) &&
					((curPos - endAp.getPos()) < windSize))
			{
				done = true;
			}
			else
			{
				window.remove();
			}
		}
		
		// And process this window
		boolean passes = this.streamWindow(streamPos, window); 
		
		// Add ourself to the window
		window.add(streamPos.currentAp);
		
		return passes;
	}

	
	/**
	 * @param streamPos
	 * @param apWind All the other APs in the window (does NOT contain the currentAp)
	 * @return
	 */
	public abstract boolean streamWindow(AlignmentPosStreamerPosition streamPos, Queue<AlignmentPos> apWind);
	
	

}
