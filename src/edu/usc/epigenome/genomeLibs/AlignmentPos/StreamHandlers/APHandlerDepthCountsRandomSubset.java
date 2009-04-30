/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;
import edu.usc.epigenome.genomeLibs.Counters.GenomicRangeCounter;
import edu.usc.epigenome.genomeLibs.Counters.TreeMapCounter;
import edu.usc.epigenome.genomeLibs.GenomicRange.*;


/**
 * @author benb
 * 
 * By default, it is a TreeMapCounter object which will randomize itself with numRandomElements upon
 * finalize.  But it can subsequently be re-used with a different number of elements by calling getDepthCounter
 * as many times as desired (it maintains a list that includes all read positions).
 *
 */
public class APHandlerDepthCountsRandomSubset extends APHandlerDepthCounts implements AlignmentPosStreamHandler {

	private static final long serialVersionUID = 4425072299373876237L;

	protected ArrayList<GenomicRange> rangeList = null;
	protected int initialArraySize = 10^7;
	
	protected int numRandomElements = 0;
	protected int numTrials = 1;
	
	/**
	 * Constructor
	 */
	public APHandlerDepthCountsRandomSubset() {
		super();
		numRandomElements = 0;
		numTrials = 1;
	}

	public APHandlerDepthCountsRandomSubset(int inNumTrials, int inNumRandomElements) {
		super();
		numRandomElements = inNumRandomElements;
		numTrials = inNumTrials;
	}

	public APHandlerDepthCountsRandomSubset(int inNumTrials, int inNumRandomElements, int inInitialArraySize) {
		super();
		this.initialArraySize = inInitialArraySize;
		numRandomElements = inNumRandomElements;
		numTrials = inNumTrials;
	}

	
	/*
	 * Access functions
	 */
	public TreeMapCounter<Integer> getDepthCounter(int maxElements)
	{
		GenomicRangeCounter ranges = getRangeCounter(maxElements);
		return ranges.getDepthCounter();
	}
	
	public GenomicRangeCounter getRangeCounter(int maxElements)
	{
		List<GenomicRange> subList = null;
				
		if (this.rangeList.size() <= maxElements)
		{
			System.err.println("Trying to take a random sublist of " + maxElements + 
					" elements, but List<GenomicRange> only has " + rangeList.size() + 
					" (will use entire list)");
			subList = rangeList;
		}
		else
		{
			// Shuffle the list
//			System.err.println("About to shuffle list of " + rangeList.size() + " elements...");
			Collections.shuffle(rangeList);
//			System.err.println("\tdone.");
			
			// Take the first maxElements elements of the shuffled list
			subList = rangeList.subList(0, maxElements);
		}
		
		// Now transform to a counter
		GenomicRangeCounter out = listToCounter(subList);
		
		return out;
	}
	
	protected static GenomicRangeCounter listToCounter(List<GenomicRange> list)
	{
		GenomicRangeCounter counter = new GenomicRangeCounter();
		counter.increment(list);
		return counter;
	}
	
	
	
	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	public void init() {
		this.rangeList = new ArrayList<GenomicRange>(initialArraySize);
	}

	public void finish() {
		
		ArrayList<TreeMapCounter<Integer>> counters = new ArrayList<TreeMapCounter<Integer>>(numTrials); 
		for (int i = 0; i < numTrials; i++)
		{
			TreeMapCounter<Integer> counter = this.getDepthCounter(numRandomElements);
			// System.err.println("Counter #" + (i+1) + "\n" + counter.excelOutput());
			counters.add(counter);
		}
		
		// Add the mean of all the randomized counters
		this.addCounts(this.meanCounter(counters));
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	public boolean streamElement(AlignmentPosStreamerPosition streamPos) 
	{
		List<GenomicRange> ranges = streamPos.currentAp.getGenomicRanges();
		this.rangeList.addAll(ranges);
		return true;
	}

	
	

}
