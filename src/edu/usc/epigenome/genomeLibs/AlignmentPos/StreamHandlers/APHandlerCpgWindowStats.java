/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

import edu.usc.epigenome.genomeLibs.CpgPair;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

/**
 * @author benb
 * 
 *
 */
/**
 * @author benb
 *
 */
abstract public class APHandlerCpgWindowStats extends APHandlerCpgHandler {

	public LinkedList<CpgPair> window;
	public int windSize = 0;
	
	public APHandlerCpgWindowStats(int inWindSize) 
	{
		windSize = inWindSize;
	}


	public void init()
	{
		window = new LinkedList<CpgPair>();
	}


	public void finish()
	{
		window = null;
	}


	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerCpgHandler#streamCpgPair(edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition, edu.usc.epigenome.genomeLibs.CpgPair)
	 * 
	 * The window does not include the CpG itself. 
	 */
	public boolean streamCpgPair(AlignmentPosStreamerPosition streamPos, CpgPair pair)
	{
		// Eat APs off the head until we are within range.
		String curChr = streamPos.currentAp.getChr();
		int curPos = streamPos.currentAp.getPos();
		boolean done = false;
		CpgPair endPair;
		while (!done && ((endPair = window.peek()) != null))
		{
			if ((endPair.getChr().equals(curChr)) &&
					((curPos - endPair.getPos()) < windSize))
			{
				done = true;
			}
			else
			{
				window.remove();
			}
		}
		
		// And process this window
		boolean passes = this.streamWindow(streamPos, pair, window); 
		
		// Add ourself to the window
		window.add(pair);
		
		return passes;
	}

	
	/**
	 * @param streamPos
	 * @param apWind All the other APs in the window (does NOT contain the currentAp)
	 * @return
	 */
	abstract public boolean streamWindow(AlignmentPosStreamerPosition streamPos, CpgPair pair, Queue<CpgPair> apWind);
	
	

}
