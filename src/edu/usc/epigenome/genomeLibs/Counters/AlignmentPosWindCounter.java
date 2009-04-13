/**
 * 
 */
package edu.usc.epigenome.genomeLibs.Counters;

import java.io.PrintStream;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosOptions;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;


/**
 * @author benb
 * 
 *
 */
public class AlignmentPosWindCounter extends GenomicRangeCounter  {


	private static final long serialVersionUID = -8312065178497620972L;

	private int windSize;
	private boolean strandSpecific;
	private AlignmentPosOptions apOptions = null;
	
	/**
	 * 
	 */
	public AlignmentPosWindCounter(int inWindSize, boolean inStrandSpecific) {
		windSize = inWindSize;
		strandSpecific = inStrandSpecific;
	}

	
	public void increment(AlignmentPos ap) {
		
		if (apOptions==null) apOptions = ap.getApOptions();
		
		if (strandSpecific)
		{
			GenomicRange fwWind = GenomicRange.generateWindowFromCoord(ap.getChr(), ap.getPos(), this.windSize, StrandedFeature.POSITIVE);
			super.increment(fwWind, ap.getDepth(true));
			GenomicRange revWind = GenomicRange.generateWindowFromCoord(ap.getChr(), ap.getPos(), this.windSize, StrandedFeature.NEGATIVE);
			super.increment(revWind, ap.getDepth(false));
		}
		else
		{
			GenomicRange wind = GenomicRange.generateWindowFromCoord(ap.getChr(), ap.getPos(), this.windSize);
			super.increment(wind, ap.getTotalDepth());
		}
	}

	

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.Counters.TreeMapCounter#excelOutput(java.lang.String, java.io.PrintStream)
	 */
	@Override
	public void excelOutput(String firstCol, PrintStream ps) {
		super.excelOutput(firstCol + "," + Integer.toString(windSize) + "," + strandSpecific + "," + apOptions.maxIdentical, ps);
	}
	
	
}
