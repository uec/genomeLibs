/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.io.PrintWriter;
import java.util.TreeMap;

import org.biojava.bio.program.gff.GFFRecord;

import edu.usc.epigenome.genomeLibs.GFFUtils;
import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;
import edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresArray;
import edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresFast;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;

/**
 * @author benb
 *
 */
public class APHandlerFeatAlignerEachfeat extends APHandlerFeatAligner {

	ChromScoresArray scoreArray = null;
	int downsamplingFactor = 1;
		
	/**
	 * @param inGtfFilename
	 * @param inWindSize
	 * @param inCensoring
	 */
	public APHandlerFeatAlignerEachfeat(String inGtfFilename, int inWindSize,
			boolean inCensoring, int inFragSize, int inDownsamplingFactor) {
		super(inGtfFilename, inWindSize, inCensoring, inFragSize);
		downsamplingFactor = inDownsamplingFactor;
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerFeatAligner#init()
	 */
	@Override
	public void init() {
		super.init();
		scoreArray = new ChromScoresArray(ChromScoresArray.ARBITRARY_GENOME, downsamplingFactor);
		scoreArray.setArbitraryGenomeLength(this.windSize*2);
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerFeatAligner#finish()
	 */
	@Override
	public void finish() {
		ListUtils.setDelim(",");
		
		
		try {
			//ChromScoresFast smoothed = scoreArray.smooth(1000, 0, 0, 100); 
			
			PrintWriter pr = new PrintWriter(System.out);
			scoreArray.singleLinePerChrOutput(pr);
			pr.flush();
		}
		catch (Exception e)
		{
			System.err.println("Could not write features:");
			e.printStackTrace(System.err);
		}
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerFeatAligner#increment(edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition, org.biojava.bio.program.gff.GFFRecord, int, double, double, double)
	 */
	@Override
	protected void increment(AlignmentPosStreamerPosition streamPos,
			String recString, int arrInd, double score, double fwScore,
			double revScore) {
		super.increment(streamPos, recString, arrInd, score, fwScore, revScore);
		
		scoreArray.addScore(recString, arrInd, new Double(score));
		//System.err.println("Adding " + score + "\t" + recString + ":" + arrInd);
	}

//	/* (non-Javadoc)
//	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerFeatStreamer#streamElement(edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition)
//	 */
//	@Override
//	public boolean streamElement(AlignmentPosStreamerPosition streamPos) {
//		return super.streamElement(streamPos);
//	}
//
//	/* (non-Javadoc)
//	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerFeatAligner#streamFeat(edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition, org.biojava.bio.program.gff.GFFRecord, int)
//	 */
//	@Override
//	protected boolean streamFeat(AlignmentPosStreamerPosition streamPos,
//			GFFRecord rec, int apRelativeOffset) {
//		return super.streamFeat(streamPos, rec, apRelativeOffset);
//	}
//
//	
	
	
}
