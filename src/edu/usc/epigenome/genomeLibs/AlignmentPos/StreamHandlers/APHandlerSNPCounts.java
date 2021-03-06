/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosSnps;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;
import edu.usc.epigenome.genomeLibs.Counters.StringCounter;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;

/**
 * @author benb
 * 
 *
 */
public class APHandlerSNPCounts extends StringCounter implements AlignmentPosStreamHandler {

	/**
	 * Constructor
	 */
	public APHandlerSNPCounts() {
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
	public boolean streamElement(AlignmentPosStreamerPosition streamPos) 
	{
		boolean passes = true;
		
		AlignmentPosSnps currentApSnps = null;
		try
		{
			currentApSnps = (AlignmentPosSnps)streamPos.currentAp;
		}
		catch (ClassCastException e)
		{
			System.err.println("APHandlerBaseCounts called with non-SNP AlignmentPos objects");
			e.printStackTrace();
			System.exit(0);
		}
		
		Iterator<ReadPos> rpIt = currentApSnps.getReadPositions().iterator();
		while (rpIt.hasNext())
		{
			// ReadPos has equals re-implemented so that identical ones are equal
			ReadPos rp = rpIt.next();
			String key = AlignmentPos.getRefTokens(streamPos.priorAps) +
			"/" + 
			streamPos.currentAp.getRefToken(true) + 
			"/" + 
			AlignmentPos.getRefTokens(streamPos.nextAps) + 
			" -> " +
			rp.getSymToken() + rp.getStrandChar();
			
//			System.err.println("\t" + key);
			
			this.increment(key);
		}
		
		return passes;
	}

	
	
	

}
