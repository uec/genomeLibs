/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

/**
 * @author benb
 * 
 * Important!  This handles both the C at postion t on the forward strand, and the C at position t+1 on the reverse strand.
 * DO NOT use with a watson-then-crick streamer.
 *
 */
 abstract public class APHandlerCpgHandler implements AlignmentPosStreamHandler {

	 public int fCpgsStreamed = 0;
	 public int fCpgsPassed = 0;
	 public int fCpgReadsStreamed = 0;
	 public int fCpgReadsPassed = 0;
	 protected boolean fReportCounts = true;
	 
	/**
	 * 
	 */
	public APHandlerCpgHandler() {
	}

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	public void init() {
	}

	public void finish() {

		if (fReportCounts)
		{
		System.err.println("Num unique CpGs streamed: " + fCpgsStreamed);
		System.err.println("Num unique CpGs streamed, passed: " + fCpgsPassed);
		System.err.println("Num CpG reads streamed: " + fCpgReadsStreamed);
		System.err.println("Num CpG reads streamed, passed: " + fCpgReadsPassed);
		}
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	public boolean streamElement(AlignmentPosStreamerPosition streamPos) 
	{
		AlignmentPosSnpsBisulfiteConverted currentBs = null;
		AlignmentPosSnpsBisulfiteConverted nextBs = null;
		try
		{
			currentBs = (AlignmentPosSnpsBisulfiteConverted)streamPos.currentAp;
			
			if ((streamPos.nextAps==null) || (streamPos.nextAps.length<1))
			{
				throw new Exception("APHandlerCpgHandler: Can not be a CpG because nextAps is empty");
			}
			else
			{
				nextBs = (AlignmentPosSnpsBisulfiteConverted)streamPos.nextAps[0].flipped();
			}
			
//			if (currentBs.getRef() != DNATools.c()) throw new Exception("APHandlerCpgHandler: First AP is not a cytosine");
//			if (nextBs.getRef() != DNATools.c()) throw new Exception("APHandlerCpgHandler: Second AP is not a cytosine");
		}
		catch (Exception e)
		{
			System.err.println("APHandlerCpgEmitter called with non-AlignmentPosSnpsBisulfiteConverted objects.");
			System.err.println("Try setting AlignmentOptions.trackBisulfiteConversion to true in AlignmentPos reader");
			e.printStackTrace();
			System.exit(0);
		}

		CpgPair pair = new CpgPair(currentBs, nextBs);
		boolean passed = streamCpgPair(streamPos, pair);
		
		// Update counters
		int depth = pair.getDepth(true) + pair.getDepth(false);
		
		 fCpgsStreamed++;
		 fCpgReadsStreamed += depth;
		 if (passed)
		 {
			 fCpgsPassed++;
			 fCpgReadsPassed += depth;
		 }
		
		
		return passed;
	}

	abstract public boolean streamCpgPair(AlignmentPosStreamerPosition streamPos, CpgPair pair); 
	

}
