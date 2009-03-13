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
 public class APHandlerCpgFeatAligner extends APHandlerCpgFeatStreamer {
	 
	 public int[] totals = null;
	 public double[] totalMeths = null;
	 
	 public double meth = -1.0;
	 public int depth = -1;

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	/**
	 * @param inGtfFilename
	 * @param inWindSize
	 */
	public APHandlerCpgFeatAligner(String inGtfFilename, int inWindSize) {
		super(inGtfFilename, inWindSize);
		// TODO Auto-generated constructor stub
	}

	
	public void init() {
		super.init();
		
		totals = new int[(this.windSize*2)+1];
		totalMeths = new double[(this.windSize*2)+1];

	}

	public void finish() {
		super.finish();

		ListUtils.setDelim(",");
		System.out.println(ListUtils.excelLine(totals));
		System.out.println(ListUtils.excelLine(totalMeths));
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


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerCpgFeatStreamer#streamFeat(edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition, edu.usc.epigenome.genomeLibs.CpgPair, org.biojava.bio.program.gff.GFFRecord, int)
	 * 
	 * Returns true if the feature uses the CpG
	 */
	protected boolean streamFeat(AlignmentPosStreamerPosition streamPos, CpgPair pair, GFFRecord rec, int cpgRelativeOffset)
	{
		AlignmentPos cur = streamPos.currentAp;
		int curPos = pair.getPos();

		boolean out = false;
		int arrInd = cpgRelativeOffset + this.windSize;

		// It's not necessarily in range
		if ((arrInd > 0) && (arrInd < totals.length)) 
		{
			// Get the meth
			if (meth == -1.0) // Cache it after the first time.
			{
//				System.err.println("\t\tGrabbing meth for " + curPos);
				meth = pair.getMethylatedFrac();
				depth = pair.getDepth(true) + pair.getDepth(false);
			}


			// *** Do we want to weight it by unique CpG or CpG reading??
			boolean COUNT_EACH_CPG_ONCE = true;
			totals[arrInd] += (COUNT_EACH_CPG_ONCE) ? 1 : depth; 
			totalMeths[arrInd] += (COUNT_EACH_CPG_ONCE) ? meth : (meth*depth);
			
			out = true;
		}

		return out;
	}

	

}
