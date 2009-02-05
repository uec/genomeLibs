/**
 * 
 */
package edu.usc.epigenome.genomeLibs.Counters;

import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;


/**
 * @author benb
 * 
 *
 */
public class GenomicWindCounter extends GenomicRangeCounter  {


	private static final long serialVersionUID = -8312065178497620972L;

	private int windSize;
	
	/**
	 * 
	 */
	public GenomicWindCounter(int inWindSize) {
		windSize = inWindSize;
	}

	public int getCount(String inChrom, int inCoord) {
		GenomicRange key = GenomicRange.generateWindowFromCoord(inChrom, inCoord, this.windSize);
		return super.getCount(key);
	}

	public void increment(String inChrom, int inCoord, int numToAdd) {
		GenomicRange key = GenomicRange.generateWindowFromCoord(inChrom, inCoord, this.windSize);
		super.increment(key, numToAdd);
	}

	public void increment(String inChrom, int inCoord) {
		GenomicRange key = GenomicRange.generateWindowFromCoord(inChrom, inCoord, this.windSize);
		super.increment(key);
	}

	public void setCount(String inChrom, int inCoord, int newCount) {
		GenomicRange key = GenomicRange.generateWindowFromCoord(inChrom, inCoord, this.windSize);
		super.setCount(key, newCount);
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
