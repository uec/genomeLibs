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
	 
	 public double score = -1.0;
	 public int depth = -1;
	 
	 public boolean doCensoring = false;

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	/**
	 * @param inGtfFilename
	 * @param inWindSize
	 */
	public APHandlerFeatAligner(String inGtfFilename, int inWindSize, boolean inCensoring) {
		super(inGtfFilename, inWindSize);
		doCensoring = inCensoring;
	}

	
	public void init() {
		super.init();
		
		totals = new int[(this.windSize*2)+1];
		totalScores = new double[(this.windSize*2)+1];
		totalScoresFw = new double[(this.windSize*2)+1];
		totalScoresRev = new double[(this.windSize*2)+1];

	}

	public void finish() {
		super.finish();

		ListUtils.setDelim(",");
		System.out.println(ListUtils.excelLine(totals));
		System.out.println(ListUtils.excelLine(totalScores));
		System.out.println(ListUtils.excelLine(totalScoresFw));
		System.out.println(ListUtils.excelLine(totalScoresRev));
	}


	


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
			 score = scoredPos.getSummaryScore();  // In the case of a CpG , this will return CpG score
			//score = streamPos.getAvgScore(true);
			depth = cur.getDepth(true) + cur.getDepth(false);

			// Strands should be relative to rec. (offset has already been flipped)
			//				System.err.println("FW Strand = " + rec.getStrand() + "\tREV strand = " + rec.getStrand().flip());
			fwScore = scoredPos.getStrandedScore(rec.getStrand());
			revScore = scoredPos.getStrandedScore(rec.getStrand().flip());

			// Now update our counters
			this.increment(streamPos, rec, arrInd, score, fwScore, revScore);
			
			out = true;
		}

		return out;
	}

	protected void increment(AlignmentPosStreamerPosition streamPos, GFFRecord rec, int arrInd, double score, double fwScore, double revScore)
	{
		totals[arrInd] += 1;
		totalScores[arrInd] += score;
		totalScoresFw[arrInd] += fwScore;
		totalScoresRev[arrInd] += revScore;
	}
	

}
