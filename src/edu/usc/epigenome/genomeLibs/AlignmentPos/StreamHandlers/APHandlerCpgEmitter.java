/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.io.File;
import java.util.*;

import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;
import edu.usc.epigenome.genomeLibs.TrackFiles.TrackFile;
import edu.usc.epigenome.genomeLibs.TrackFiles.TrackFileRandomAccess;

/**
 * @author benb
 * 
 *
 */
 public class APHandlerCpgEmitter extends APHandlerCpgHandler {
	 
	 public String cpgTrackFilename = null;
	 public TrackFile cpgTrackFile = null;
	 public int cpgDensWindsize = 200;

	/**
	 * 
	 */
	public APHandlerCpgEmitter() {
	}

	public APHandlerCpgEmitter(String inCpgDensityTrackFilename) {
		cpgTrackFilename = inCpgDensityTrackFilename;
	}

	public APHandlerCpgEmitter(String inCpgDensityTrackFilename, int inWindSize) {
		cpgTrackFilename = inCpgDensityTrackFilename;
		cpgDensWindsize = inWindSize;
	}

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	public void init() {

		// Open our trackfile
		if (cpgTrackFilename!=null)
		{
			try
			{
				cpgTrackFile = new TrackFileRandomAccess(new File(cpgTrackFilename), "hg18", false, false);
			}
			catch (Exception e)
			{
				System.err.println("Can't open CpG track file " + this.cpgTrackFilename + ":\n" + e.toString());
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		System.out.println(headers());
		super.init();
	}

	public void finish() {
		
		// Close our trackfile
		
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

	public boolean streamCpgPair(AlignmentPosStreamerPosition streamPos, CpgPair pair)
	{
		// PRINT THE MAXIDENTICAL DEPTH AND THE NON 
		
		AlignmentPos cur = streamPos.currentAp;
		
		// Put a number so matlab can open it.
		long globalPos = GoldAssembly.getGlobalOffset(cur.getChr(), "hg18", cur.getPos());
//		System.out.print(cur.getChr() + ",");
		System.out.print(globalPos + ",");
		System.out.print(cur.getPos() + ",");
		
		
//		System.out.print( streamPos.priorAps.length + ",");
////		System.out.print( streamPos.preNmerCpgDensity() + ",");
//		System.out.print( AlignmentPos.getCpgDensityStr(streamPos.priorAps) + ",");
//		System.out.print( streamPos.nextAps.length + ",");
////		System.out.print( streamPos.nextNmerCpgDensity() + ",");
//		System.out.print( AlignmentPos.getCpgDensityStr(streamPos.nextAps) + ",");
		
		double preDens = -1.0, postDens = -1.0;
		if (cpgTrackFile!=null)
		{
			try
			{
				preDens = cpgTrackFile.getValAvg(cur.getChr(), cur.getPos()-cpgDensWindsize, cpgDensWindsize);
				postDens = cpgTrackFile.getValAvg(cur.getChr(), cur.getPos(), cpgDensWindsize);
			}
			catch (Exception e)
			{
				System.err.println("Can't read CpG track file " + this.cpgTrackFilename + ":\n" + e.toString());
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		System.out.print( cpgDensWindsize + ",");
		System.out.print( preDens + ",");
		System.out.print( cpgDensWindsize + ",");
		System.out.print( postDens + ",");
		
		// Add ourself to the window since we aren't included
		System.out.print(pair.csvStats());
		System.out.println();
		
		return true;
	}
	

}
