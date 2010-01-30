package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.util.LinkedList;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgSummarizer;

/**
 * @author benb
 *
 * Walks over a sliding window of adjacent Cpgs
 */

public abstract class CpgWalker {


/*
 * 
 * 	
 */
	public CpgWalkerParams walkParams = null;
	
	// List management
	public LinkedList<Cpg> window = new LinkedList<Cpg>();
	
	// Some useful summarizers for the window
	protected CpgSummarizer methSummarizer = new CpgMethLevelSummarizer();

	
	
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
		methSummarizer.streamCpg(cpg);
		
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
				methSummarizer.removeCpg(endCpg);
			}
		}
		
		// And process the window
		if (window.size()>walkParams.minCpgs)
		{
			this.processWindow();
		}
		
	}

	// OVERRIDE THESE
	abstract protected void processWindow();

	
}
