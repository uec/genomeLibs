/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.Counters.StringCounter;

/**
 * @author benb
 * 
 *
 */
public class APHandlerGffEmitter extends StringCounter implements AlignmentPosStreamHandler {

	/**
	 * 
	 */
	public APHandlerGffEmitter() {
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
	public boolean streamElement(AlignmentPos[] priorAps,
			AlignmentPos currentAp, AlignmentPos[] nextAps) 
	{
		boolean passes = true;
		System.out.print(currentAp.gffLine());
		return passes;
	}

	
	
	

}
