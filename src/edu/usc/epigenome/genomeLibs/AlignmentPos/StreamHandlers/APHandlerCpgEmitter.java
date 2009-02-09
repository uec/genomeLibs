/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.*;

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



	public boolean streamCpgPair(AlignmentPos[] priorAps,
			CpgPair pair, AlignmentPos[] nextAps)
	{
		System.out.print( priorAps.length + ",");
		System.out.print( AlignmentPos.getCpgDensityStr(priorAps) + ",");
		System.out.print( nextAps.length + ",");
		System.out.print( AlignmentPos.getCpgDensityStr(nextAps) + ",");
		System.out.print(pair.csvLine());
		System.out.println();
		
		return true;
	}
	

}
