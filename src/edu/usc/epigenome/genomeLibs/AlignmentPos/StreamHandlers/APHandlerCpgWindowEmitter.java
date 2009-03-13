/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import BisulfiteCytosines.CpgPair;

import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

/**
 * @author benb
 * 
 *
 */
 public class APHandlerCpgWindowEmitter extends APHandlerCpgWindowStats {


	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	/**
	 * @param inWindSize
	 */
	public APHandlerCpgWindowEmitter(int inWindSize) {
		super(inWindSize);
	}

	public void init() {
		System.out.println(headers());
		super.init();
	}

	public void finish() {
		super.finish();
	}

	public static String headers()
	{
		String out = "";
		
		out += "chr,";
		out += "pos,";
		out += "preWindLen,";
		out += "preWindCpgDens,";
		out += "postWindLen,";
		out += "postWindCpgDens,";
		out += CpgPair.csvStatsHeaders();
		
		return out;
	}

	public boolean streamWindow(AlignmentPosStreamerPosition streamPos, CpgPair pair, Queue<CpgPair> cpgWind)
	{

		System.out.print(streamPos.currentAp.getChr() + ",");
		System.out.print(streamPos.currentAp.getPos() + ",");
		
		System.out.print( streamPos.priorAps.length + ",");
//		System.out.print( streamPos.preNmerCpgDensity() + ",");
		System.out.print( AlignmentPos.getCpgDensityStr(streamPos.priorAps) + ",");
		System.out.print( streamPos.nextAps.length + ",");
//		System.out.print( streamPos.nextNmerCpgDensity() + ",");
		System.out.print( AlignmentPos.getCpgDensityStr(streamPos.nextAps) + ",");
		
		// Add ourself to the window since we aren't included
		Vector<CpgPair> list = new Vector<CpgPair>(cpgWind);
		list.add(pair);
		System.out.print(CpgPair.csvStats(list));
		System.out.println();
		
		return true;
	}


}
