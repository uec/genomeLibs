/**
 * 
 */
package edu.usc.epigenome.genomeLibs.Counters;

import java.io.PrintStream;
import java.util.Iterator;

import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;


/**
 * @author benb
 * 
 * Very similar to GenomicRangeCounter, except that this one expands ranges to 
 * the genomic "windows" they fall into , which are defined as non-overlapping
 * segments of windSize across the genome.
 * 
 * ***** TO DO .  Right now i just use the start coordinate of the input 
 * ***** range, whereas we should probably increment all overlapping windows.
 */
public class GenomicRangeWindCounter extends GenomicRangeCounter  {

	private static final long serialVersionUID = 527845714730797504L;


	/* Object methods */
	protected int windSize;
	protected boolean strandSpecific;
	protected String genome;
	
	/**
	 * 
	 */
	public GenomicRangeWindCounter(int inWindSize, boolean inStrandSpecific, String inGenome) {
		windSize = inWindSize;
		strandSpecific = inStrandSpecific;
		genome = inGenome;
	}

	
	/* Readers */
	public int getCount(GenomicRange inRange) {
		return super.getCount(windFromRange(inRange));
	}

	
	
	/* Writers */
	public void increment(GenomicRange inRange) {
		super.increment(windFromRange(inRange));
	}

	public void increment(GenomicRange inRange, int numToAdd) {
		super.increment(windFromRange(inRange), numToAdd);
	}

	public void setCount(GenomicRange inRange, int newCount) {
		super.setCount(windFromRange(inRange), newCount);
	}
	
	/* Utility functions */

	protected GenomicRange windFromRange(GenomicRange inRange)
	{
		return (this.strandSpecific) ? 
				GenomicRange.generateWindowFromCoord(inRange.getChrom(), inRange.getStart(), this.windSize, inRange.getStrand()) : 
					GenomicRange.generateWindowFromCoord(inRange.getChrom(), inRange.getStart(), this.windSize);
	}
	
	/* Output */
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.Counters.TreeMapCounter#excelOutput(java.lang.String, java.io.PrintStream)
	 */
	@Override
	public void excelOutput(String firstCol, PrintStream ps) {

		String initialCols = firstCol + "," + Integer.toString(windSize) + "," + strandSpecific;

		// Output all windows
		Iterator<GenomicRange> rangeIt = null;
		try {
			rangeIt = GenomicRange.allPossibleGenomicRanges(genome, windSize, strandSpecific, true);
		} catch (Exception e) {
			System.err.println("GenomicRange.allPossibleGenomicRanges threw an exception: ");
			e.printStackTrace();
			System.exit(0);
		}
		
		while (rangeIt.hasNext())
		{
			GenomicRange key = rangeIt.next();
			int count = super.getCount(key);

			ps.print(initialCols);
			ps.print("," + key.toString());
			ps.print("," + count);
			ps.print("\n");
		}
	}
	
	
}
