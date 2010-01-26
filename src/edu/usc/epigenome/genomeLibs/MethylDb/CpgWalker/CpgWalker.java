package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.util.LinkedList;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

/**
 * @author benb
 *
 * Walks over a sliding window of adjacent Cpgs
 */

public abstract class CpgWalker {


/*
 * 
 * Uses a CpgWalkerParams class to store parameters
 * stepSize: We do a simple wiggle track with values every n bp (50 or 100)
 * windSize: 
 * minCpgs: The window averages this many Cpgs, centered on the step position.
 *   If there are more than minCpgs withing the windSize, we use all of them.
 * maxWindStretch: We only stretch the window this far in either direction to
 *   look for Cpgs to average.  If there are none, we put a NaN value.
 *   
 * Do coverage as well as meth
 * Do deltas in coverage and meth.
 * 	
 */
	public CpgWalkerParams walkParams = null;
	
	// List management
	public LinkedList<Cpg> window = new LinkedList<Cpg>();

	
	
	/**
	 * 
	 */
	public CpgWalker(CpgWalkerParams inWalkParams) {
		super();
		this.walkParams = inWalkParams;
	}
	
	public void reset()
	{
		window = new LinkedList<Cpg>();
	}
	
	/**
	 * @param cpg
	 */
	public void streamCpg(Cpg cpg)
	{
		int newPos = cpg.chromPos;
		
		// Add this Cpg to the head of the queue
		window.add(cpg);
		
		// Remove cpgs from the tail
		boolean done = false;
		Cpg endCpg;
		while (!done && ((endCpg = window.peek()) != null))
		{
			if ((newPos - endCpg.chromPos) < this.walkParams.maxWindSize)
			{
				done = true;
			}
			else
			{
				window.remove();
			}
		}
		
		// And process the window
		if (window.size()>0)
		{
			this.processWindow();
		}
		
	}

	// OVERRIDE THESE
	abstract protected void processWindow();

	
}
