package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.io.PrintWriter;

import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;

/**
 * @author benb
 *
 * Walks over a sliding window of adjacent Cpgs
 */

public class CpgWalker {


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
	protected PrintWriter outStream = null;
	
	// List management
	
	
	
	/**
	 * 
	 */
	public CpgWalker(CpgWalkerParams inWalkParams,
			PrintWriter outStream) {
		super();
		this.walkParams = inWalkParams;
	}

	
}
