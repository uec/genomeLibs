/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;
import edu.usc.epigenome.genomeLibs.Counters.ReadPosCounter;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;

/**
 * @author benb
 * 
 *
 */
public class APHandlerSymbolCounts extends ReadPosCounter implements AlignmentPosStreamHandler {

	/**
	 * 
	 */
	private static final long serialVersionUID = -3266222243481309526L;



	/**
	 *  Constructor
	 */
	public APHandlerSymbolCounts() {
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
		Iterator<ReadPos> it = streamPos.currentAp.getReadPositions().iterator();
		while (it.hasNext())
		{
			ReadPos rp = it.next();
			this.increment(rp);
		}
		return passes;
	}



	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TreeMapCounter#excelOutput()
	 */
	@Override
	public String excelOutput() {
		// TODO Auto-generated method stub
		return super.excelOutput();
	}
}
