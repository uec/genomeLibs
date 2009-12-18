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
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicPositionScored;
import edu.usc.epigenome.genomeLibs.TrackFiles.TrackFile;
import edu.usc.epigenome.genomeLibs.TrackFiles.TrackFileRandomAccess;

/**
 * @author benb
 * 
 *
 */
 public class APHandlerFeatAligner extends APHandlerFeatStreamer {
	 
	 public int[] totals = null;
	 public double[] totalScores = null;
	 public double[] totalScoresFw = null;
	 public double[] totalScoresRev = null;
	 public int arrLen;
	 
	 public double score = -1.0;
	 public int depth = -1;
	 
	 public boolean doCensoring = false;
	 public int fragSize = 500;

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	/**
	 * @param inGtfFilename
	 * @param inWindSize
	 */
	public APHandlerFeatAligner(String inGtfFilename, int inWindSize, boolean inCensoring, int inFragSize) {
		super(inGtfFilename, inWindSize);
		doCensoring = inCensoring;
		fragSize = inFragSize;
	}

	
	public void init() {
		super.init();
		
		totals = new int[(this.windSize*2)+1];
		totalScores = new double[(this.windSize*2)+1];
		totalScoresFw = new double[(this.windSize*2)+1];
		totalScoresRev = new double[(this.windSize*2)+1];
		arrLen = totals.length;
	}

	public void finish() {
		super.finish();

		ListUtils.setDelim(",");
		System.out.println(ListUtils.excelLine(totals));
		System.out.println(ListUtils.excelLine(totalScores));
		System.out.println(ListUtils.excelLine(totalScoresFw));
		System.out.println(ListUtils.excelLine(totalScoresRev));
	}


	


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerFeatStreamer#streamFeat(edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition, org.biojava.bio.program.gff.GFFRecord, int)
	 * 
	 * We use GFFUtils.gffCsvPosLine(rec) to generate a unique record string
	 */
	@Override
	protected boolean streamFeat(AlignmentPosStreamerPosition streamPos, GFFRecord rec, int apRelativeOffset)
	{
		AlignmentPos cur = streamPos.currentAp;
		GenomicPositionScored scoredPos = streamPos.currentScoredPosition;
		int curPos = cur.getPos();

		boolean out = false;
		int arrInd = apRelativeOffset + this.windSize;

		// It's not necessarily in range
		if ((arrInd > 0) && (arrInd < totals.length)) 
		{
			// Get the score
			double fwScore = -1.0;
			double revScore = -1.0;

			// Get the summary scores
			// score = scoredPos.getSummaryScore();  // In the case of a CpG , this will return CpG score
			score = streamPos.getAvgScore(true); //rec.getStrand());
			depth = cur.getDepth(true) + cur.getDepth(false);

			// Strands should be relative to rec. (offset has already been flipped)
			//				System.err.println("FW Strand = " + rec.getStrand() + "\tREV strand = " + rec.getStrand().flip());
			fwScore = scoredPos.getStrandedScore(rec.getStrand());
			revScore = scoredPos.getStrandedScore(rec.getStrand().flip());

			// Make a string for the rec
			//String recString = GFFUtils.gffCsvPosLine(rec) + "," + GFFUtils.getGffRecordName(rec);
			String recString = GFFUtils.gffCsvPosLine(rec);
			
			// Now update our counters
			if (cur.getApOptions().onlyFirstCycle)
			{

				// EVEN THOUGH IT'S MORE EFFICIENT, THIS ONE LOOKED TOO CHOPPY
//				// We are only streaming the read start position(s).
//				// Increment the midpoint of the fragment based on fragment length.
//				if (fwScore > 0.0)
//				{
//					int fwPos = arrInd+Math.round(this.fragSize/2);
//					if (fwPos < arrLen) increment(streamPos, recString, fwPos, fwScore, fwScore, 0.0);
//					int revPos = arrInd-Math.round(this.fragSize/2);
//					if (revPos >= 0) increment(streamPos, recString, revPos, revScore, 0.0, revScore);
//				}
				
				
				
				// We are only streaming the read start position(s). So we will go through and 
				// increment every position in the read (using the fragLength param)
				int stPos, endPos, i;
				if (this.fragSize<2) this.fragSize=2;
				if (fwScore > 0.0)
				{
					stPos = arrInd;
					endPos = Math.min(arrLen-1, arrInd+(this.fragSize/2)-1);
					//System.err.println("\tIncrementing " + GFFUtils.getGffRecordName(rec) + "\t FW from " + stPos + " to " + endPos + "\t" + fwScore + "\tfrag=" + this.fragSize);
					for (i=stPos; i<=endPos; i++) this.increment(streamPos, recString, i, fwScore, fwScore, 0.0);
				}
				if (revScore > 0.0)
				{
					stPos = Math.max(0, arrInd-(this.fragSize/2)+1);
					endPos = arrInd;
					//System.err.println("\tIncrementing " + GFFUtils.getGffRecordName(rec) + "\t REV from " + stPos + " to " + endPos  + "\t" + revScore);
					for (i=stPos; i<=endPos; i++) this.increment(streamPos, recString, i, revScore, 0.0, revScore);
				}
			}
			else
			{
				// Otherwise, just increment the current position
				this.increment(streamPos, recString, arrInd, score, fwScore, revScore);
			}
			
			out = true;
		}

		return out;
	}

	// Provide a unique string for the record so that it will be easily hashable
	protected void increment(AlignmentPosStreamerPosition streamPos, String uniqueRecordString, int arrInd, double score, double fwScore, double revScore)
	{
		totals[arrInd] += 1;
		totalScores[arrInd] += score;
		totalScoresFw[arrInd] += fwScore;
		totalScoresRev[arrInd] += revScore;
	}
	

}
