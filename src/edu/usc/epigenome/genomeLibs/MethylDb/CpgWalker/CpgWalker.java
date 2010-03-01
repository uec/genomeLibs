package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.util.LinkedList;
import java.util.logging.Logger;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizerStrandspecific;
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
	protected CpgSummarizer methSummarizerFw = new CpgMethLevelSummarizerStrandspecific(true);
	protected CpgSummarizer methSummarizerRev = new CpgMethLevelSummarizerStrandspecific(false);
	public String lastChrom = "noChrom";

	
	
	/**
	 * 
	 */
	public CpgWalker(CpgWalkerParams inWalkParams) {
		super();
		this.walkParams = inWalkParams;
//		System.out.printf("track name=\"%s\" description=\"%s\" useScore=0 itemRgb=On visibility=4\n",
//				"test", "test");

	}
	
	public void newChrom()
	{
		reset();
	}
	
	public void reset()
	{
		window = new LinkedList<Cpg>();
		methSummarizer = new CpgMethLevelSummarizer();
		methSummarizerFw = new CpgMethLevelSummarizerStrandspecific(true);
		methSummarizerRev = new CpgMethLevelSummarizerStrandspecific(false);
	}
	
	/**
	 * @param cpg
	 */
	public void streamCpg(Cpg cpg)
	{
		int newPos = cpg.chromPos;
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(
				String.format("Found Cpg: %d\n", newPos));
		
		// Add this Cpg to the head of the queue
		window.add(cpg);
		methSummarizer.streamCpg(cpg);
		methSummarizerFw.streamCpg(cpg);
		methSummarizerRev.streamCpg(cpg);
		
		
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
				methSummarizerFw.removeCpg(endCpg);
				methSummarizerRev.removeCpg(endCpg);
			}
		}
		
		//System.err.println("\tChecking " + this.windStr());
		
		// And process the window
		if (window.size()>=walkParams.minCpgs)
		{
//			System.err.println("\t\t Sufficient Cpgs");
//			double mean = this.methSummarizer.getValMean(true);
//			System.err.println("\t\t Mean meth=" + mean);
//			if (mean < 0.7)
//			System.out.println(MethylDbUtils.bedLine("chr11", windStart(), windEnd(), ".", mean));
			this.processWindow();
		}
		else
		{
//			System.err.println("\t\t Not enough Cpgs");
		}
		
	}

	public String windStr()
	{
		return windStr(false);
	}
	
	public String windStr(boolean longVers)
	{
		StringBuffer sb = new StringBuffer(10000);
		
		int e = windEnd();
		int s = windStart();
		sb.append(String.format("Wind %d-%d (%d bp) has %d elements: ",
				s, e, e-s+1, this.window.size()));
		if (longVers)
		{
			for (int i = 0; i < window.size(); i++)
			{
				if (i>0) sb.append(", ");
				Cpg cpg = window.get(i);
				sb.append(cpg.chromPos);
			}
		}
		
		
		return sb.toString();
	}
	
	public int windStart()
	{
		if (window.size()==0) return -1;
		return window.get(0).chromPos;
	}

	public int windEnd()
	{
		if (window.size()==0) return -1;
		return window.get(window.size()-1).chromPos;
	}
	
	// OVERRIDE THESE
	abstract protected void processWindow();

	
}
