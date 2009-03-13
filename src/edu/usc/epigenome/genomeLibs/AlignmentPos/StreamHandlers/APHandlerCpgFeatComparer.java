/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.io.File;
import java.util.*;

import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;
import org.usckeck.genome.ChromFeatures;

import BisulfiteCytosines.CpgPair;

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
 public class APHandlerCpgFeatComparer extends APHandlerCpgFeatStreamer {
	 
	 public double meth = -1.0;
	 public int depth = -1;

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	/**
	 * @param inGtfFilename
	 * @param inWindSize
	 */
	public APHandlerCpgFeatComparer(String inGtfFilename, int inWindSize) {
		super(inGtfFilename, inWindSize);
		// TODO Auto-generated constructor stub
	}

	
	public void init() {
		super.init();
		
		// Make the header
		System.out.println(headers());
		
	}

	public void finish() {
		super.finish();
	}


	public static String headers()
	{
		String out = "";
		
		out += GFFUtils.gffCsvMinimalLineHeader();
		out += ",chr,cpgPos,";
		out += CpgPair.csvStatsHeaders();
		return out;
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerCpgFeatStreamer#streamCpgPair(edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition, edu.usc.epigenome.genomeLibs.CpgPair)
	 */
	@Override
	public boolean streamCpgPair(AlignmentPosStreamerPosition streamPos,
			CpgPair pair) {

		// Set out caches up
		meth = -1.0;
		depth = -1;
		
		return super.streamCpgPair(streamPos, pair);
	}


	protected boolean streamFeat(AlignmentPosStreamerPosition streamPos, CpgPair pair, GFFRecord rec, int cpgRelativeOffset)
	{
		AlignmentPos cur = streamPos.currentAp;
		int curPos = pair.getPos();
		String chr = cur.getChr();

		boolean out = true;

		System.out.print(GFFUtils.gffCsvMinimalLine(rec));
		System.out.print(",");
		System.out.print(chr);
		System.out.print(",");
		System.out.print(curPos);
		System.out.print(",");
		
		// Finally the csv stats
		System.out.print(pair.csvStats());
		System.out.println();

		
		return out;
	}

	

}
