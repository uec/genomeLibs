/**
 * 
 */
package edu.usc.epigenome.genomeLibs.Counters;

import java.util.Iterator;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;


/**
 * @author benb
 * 
 *
 */
public class GenomicRangeCounter extends TreeMapCounter<GenomicRange>  {


	private static final long serialVersionUID = -3437283928405335060L;

	/**
	 * 
	 */
	public GenomicRangeCounter() {
	}
	
	public TreeMapCounter<Integer> getDepthCounter()
	{
		TreeMapCounter<Integer> counts = new TreeMapCounter<Integer>();
		
		Iterator<GenomicRange> it = this.keySet().iterator();
		while (it.hasNext())
		{
			GenomicRange range = it.next();
			int depth = this.getCount(range);
			
			if (range.getStrand() == StrandedFeature.NEGATIVE) depth *= -1;
			counts.increment(new Integer(depth));
		}
		
		return counts;
	}
	
	
	
}
