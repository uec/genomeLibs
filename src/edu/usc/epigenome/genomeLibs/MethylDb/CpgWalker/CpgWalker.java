package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.TabularOutput;
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

public abstract class CpgWalker implements TabularOutput {


/*
 * 
 * 	
 */
	public CpgWalkerParams walkParams = null;
	protected boolean useSummarizers = true;
	
	// List management
	private LinkedList<Cpg> window = new LinkedList<Cpg>();
	
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
		this.resetSummarizers();
	}
	
	public void resetSummarizers()
	{
		if (useSummarizers)
		{
			methSummarizer = new CpgMethLevelSummarizer();
			methSummarizerFw = new CpgMethLevelSummarizerStrandspecific(true);
			methSummarizerRev = new CpgMethLevelSummarizerStrandspecific(false);
		}
	}
	
	/**
	 * This uses a fixed window size and is much faster
	 * @param cpg Must be streamed in serially.
	 */
	public void streamCpgFixedWind(Cpg cpg)
	{
		int newPos = cpg.chromPos;
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(
				String.format("Found Cpg: %d\n", newPos));

		// Add this Cpg to the head of the queue
		window.add(cpg);
		if (useSummarizers)
		{
			methSummarizer.streamCpg(cpg);
			methSummarizerFw.streamCpg(cpg);
			methSummarizerRev.streamCpg(cpg);
		}


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
				if (useSummarizers)
				{
					methSummarizer.removeCpg(endCpg);
					methSummarizerFw.removeCpg(endCpg);
					methSummarizerRev.removeCpg(endCpg);
				}
			}
		}
		
		//System.err.println("\tChecking " + this.windStr());
		
		// And process the window
		if (window.size()>=walkParams.minCpgs)
		{
//			System.err.println("\t\t Sufficient Cpgs: " + window.size());
//			double mean = this.methSummarizer.getValMean(true);
//			System.err.println("\t\t Mean meth=" + mean);
//			if (mean < 0.7)
//			System.out.println(MethylDbUtils.bedLine("chr11", windStart(), windEnd(), ".", mean));
			this.processWindow(this.window);
		}
		else
		{
//			System.err.println("\t\t Not enough Cpgs");
		}
	}
		
	/**
	 * This uses a fixed step, and can have variable window sizing to accommodate
	 * widely varying CG density.
	 * @param cpg Must be streamed in serially.
	 */
	public void streamCpgVariableWind(Cpg cpg)
	{
		int newPos = cpg.chromPos;
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(
				String.format("Variable wind found Cpg: %d\n", newPos));

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

		//System.err.println("\tChecking " + this.windStr());
		
		// And process the window
		//System.err.println(this.windStr());
		if (window.size()>=walkParams.minCpgs)
		{
			
			// First we set our summarizers.
			// Take the last minCpgs as the sub-window.  Do it as an iterator since
			// a linked list might be quicker iterating than using "get(i)".
			Iterator<Cpg> backIt = window.descendingIterator();
			int i = 0;
			int lastPos = 0;
			if (useSummarizers) this.resetSummarizers();
			while (i<walkParams.minCpgs)
			{
				// ***** REMOVE THIS
				if (!backIt.hasNext()) System.err.println("Why did we run out of window elements?!"); // REMOVE
				// ***** REMOVE THIS

				Cpg backCpg = backIt.next();
				//System.err.println("\t(i=" + i + ") cpg=" + backCpg.chromPos);
				lastPos = backCpg.chromPos;
				
				if (useSummarizers)
				{
					methSummarizer.streamCpg(backCpg);
					methSummarizerFw.streamCpg(backCpg);
					methSummarizerRev.streamCpg(backCpg);
				}
				
				i++;
			}
			
			
			
			
//			System.err.println("\tDomain size=" + (cpg.chromPos-lastPos+1));
//			System.out.println("chr11\t" + lastPos + "\t" + cpg.chromPos);
			
			
			// Process this window
			this.processWindow(this.window.subList(this.window.size()-walkParams.minCpgs, this.window.size()));
		}
	}

	/**
	 * @param cpg
	 */
	public void streamCpg(Cpg cpg)
	{
		//System.err.println("useFixedStep=" + walkParams.useFixedStep);
		if (this.walkParams.useVariableWindow)
		{
			this.streamCpgVariableWind(cpg);
		}
		else
		{
			this.streamCpgFixedWind(cpg);
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
	abstract protected void processWindow(List<Cpg> inWindow);

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TabularOutput#headerStr()
	 */
	@Override
	public String headerStr() throws Exception {
		return null;
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TabularOutput#toCsvStr()
	 */
	@Override
	public String toCsvStr() throws Exception {
		return null;
	}

	
}
