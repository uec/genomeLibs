/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

/**
 * @author benb
 * 
 *
 */
 public class APHandlerCpgEmitter extends APHandlerCpgHandler {

	/**
	 * 
	 */
	public APHandlerCpgEmitter() {
	}

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	public void init() {
	}

	public void finish() {
	}



	public boolean streamCpgPair(AlignmentPosStreamerPosition streamPos, CpgPair pair)
	{
		System.out.print( streamPos.priorAps.length + ",");
		System.out.print( streamPos.preNmerCpgDensity() + ",");
//		System.out.print( AlignmentPos.getCpgDensityStr(streamPos.priorAps) + ",");
		System.out.print( streamPos.nextAps.length + ",");
		System.out.print( streamPos.nextNmerCpgDensity() + ",");
//		System.out.print( AlignmentPos.getCpgDensityStr(streamPos.nextAps) + ",");
		System.out.print(pair.csvLine());
		System.out.println();
		
		return true;
	}
	

}
