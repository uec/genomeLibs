/**
 * 
 */
package edu.usc.epigenome.genomeLibs.Counters;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
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
	
	/**
	 * 
	 */
	public AlignmentPosWindCounter(int inWindSize, boolean inStrandSpecific) {
		windSize = inWindSize;
		strandSpecific = inStrandSpecific;
	}

	
	public void increment(AlignmentPos ap) {
		
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
	 * @see edu.usc.epigenome.genomeLibs.TreeMapCounter#excelOutput()
	 */
	@Override
	public String excelOutput() {
		return super.excelOutput(Integer.toString(windSize));
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TreeMapCounter#excelOutput(java.lang.String)
	 */
	@Override
	public String excelOutput(String firstCol) {
		return super.excelOutput(firstCol + "," + Integer.toString(windSize));
	}
	
	
	
}
