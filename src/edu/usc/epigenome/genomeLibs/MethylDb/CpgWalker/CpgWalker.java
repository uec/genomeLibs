package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.TabularOutput;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizerStranded;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgSummarizer;

/**
 * @author benb
 *
 * Walks over a sliding window of adjacent Cpgs
 */

public class CpgWalker implements TabularOutput {

/*
 * 
 * 	
 */
	public CpgWalkerParams walkParams = null;
	protected boolean useSummarizers = true;
	public static final int PROCESS_WINDOW_EVENT = 1;
	protected List<Cpg> lastProcessedWindow = null;
	
	// List management
	private LinkedList<Cpg> window = new LinkedList<Cpg>();
	
	// Some useful summarizers for the window
	protected CpgSummarizer methSummarizer = new CpgMethLevelSummarizer();
	protected CpgSummarizer methSummarizerFw = new CpgMethLevelSummarizerStranded(true);
	protected CpgSummarizer methSummarizerRev = new CpgMethLevelSummarizerStranded(false);
	String curChr = null;
	public String lastChrom = "noChrom";

	// Sends events when it gets a good window
	protected List<ActionListener> listeners = new ArrayList<ActionListener>(10);
	
	/**
	 * 
	 */
	public CpgWalker(CpgWalkerParams inWalkParams) {
		super();
		this.walkParams = inWalkParams;
//		System.out.printf("track name=\"%s\" description=\"%s\" useScore=0 itemRgb=On visibility=4\n",
//				"test", "test");

	}
	
	public void clearWindowListeners()
	{
		listeners = new ArrayList<ActionListener>(10);
	}
	
	public void addWindowListener(ActionListener l)
	{
		listeners.add(l);
	}
	
	/**
	 * @return the curChr
	 */
	public String getCurChr() {
		return curChr;
	}

	/**
	 * @param curChr the curChr to set
	 */
	public void setCurChr(String curChr) {
		this.curChr = curChr;
		reset();
		this.alertNewChrom();
	}

	protected boolean onNewChrom()
	{
		return !this.lastChrom.equalsIgnoreCase(this.getCurChr());
	}


	/**
	 * Outputs the last domain
	 */
	public void finishChr() {
	}

	
	protected void alertNewChrom()
	{
	}
	
	public void reset()
	{
		window = new LinkedList<Cpg>();
		this.resetSummarizers();
		System.err.println("CpgWalker reset()");
	}
	
	public void resetSummarizers()
	{
		if (useSummarizers)
		{
			methSummarizer = new CpgMethLevelSummarizer();
			methSummarizerFw = new CpgMethLevelSummarizerStranded(true);
			methSummarizerRev = new CpgMethLevelSummarizerStranded(false);
		}
	}
	
	public List<Cpg> getLastProcessedWindow() {
		return lastProcessedWindow;
	}

	/**
	 * This uses a fixed window size and is much faster
	 * @param cpg Must be streamed in serially.
	 */
	protected void streamCpgFixedWind(Cpg cpg)
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
			if ((newPos - endCpg.chromPos) < this.walkParams.maxScanningWindSize)
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
		if (window.size()>=walkParams.minScanningWindCpgs)
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
	protected void streamCpgVariableWind(Cpg cpg)
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
			if ((newPos - endCpg.chromPos) < this.walkParams.maxScanningWindSize)
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
		if (window.size()>=walkParams.minScanningWindCpgs)
		{
			
			// First we set our summarizers.
			// Take the last minCpgs as the sub-window.  Do it as an iterator since
			// a linked list might be quicker iterating than using "get(i)".
			Iterator<Cpg> backIt = window.descendingIterator();
			int i = 0;
			int lastPos = 0;
			if (useSummarizers) this.resetSummarizers();
			boolean minSizeReached = false;
			while ((i<window.size()) && !(minSizeReached && (i>walkParams.minScanningWindCpgs)))
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
				
				int windLen = newPos-lastPos;
				minSizeReached = (windLen>=this.walkParams.minScanningWindSize);
				//System.err.printf("Checking size %d (minSizeReached=%s)\n",windLen,minSizeReached);

				i++;
			}


			// Process this window
	//		System.err.printf("Checking size %d (minSizeReached=%s)\n",cpg.chromPos-lastPos+1,minSizeReached);
			if (minSizeReached && (i>walkParams.minScanningWindCpgs))
			{
				this.processWindow(this.window.subList(this.window.size()-(i-1), this.window.size()));
			}
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

	public String lastProcessedWindStr(boolean longVers)
	{
		return windStr(this.getLastProcessedWindow(),longVers);
	}

	public String windStr()
	{
		return windStr(this.window);
	}
	
	public static String windStr(List<Cpg> inWind)
	{
		return windStr(inWind, false);
	}

	public String windStr(boolean longVers)
	{
		return windStr(this.window, longVers);
	}

	
	public static String windStr(List<Cpg> inWind, boolean longVers)
	{
		StringBuffer sb = new StringBuffer(10000);
		
		int e = windEnd(inWind);
		int s = windStart(inWind);
		sb.append(String.format("Wind %d-%d (%d bp) has %d elements: ",
				s, e, e-s+1, inWind.size()));
		if (longVers)
		{
			for (int i = 0; i < inWind.size(); i++)
			{
				if (i>0) sb.append(", ");
				Cpg cpg = inWind.get(i);
				sb.append(cpg.chromPos);
			}
		}
		
		
		return sb.toString();
	}

	public int windLen()
	{
		return windLen(this.window);
	}
	
	public int windStart()
	{
		return windStart(window);
	}

	public int windEnd()
	{
		return windEnd(window);
	}
	
	public static int windLen(List<Cpg> inWind)
	{
		return windEnd(inWind) - windStart(inWind) + 1;
	}
	
	public static int windStart(List<Cpg> inWind)
	{
		if (inWind.size()==0) return -1;
		return inWind.get(0).chromPos;
	}

	public static int windEnd(List<Cpg> inWind)
	{
		if (inWind.size()==0) return -1;
		return inWind.get(inWind.size()-1).chromPos;
	}

	// !!! OVERRIDE THIS !!!
	protected void processWindow(List<Cpg> inWindow)
	{
		lastProcessedWindow = inWindow;
		
		// Basic version does nothing except alert listeners
		ActionEvent e = new ActionEvent(this, PROCESS_WINDOW_EVENT, null);
		for (ActionListener l : listeners)
		{
			l.actionPerformed(e);
		}

		this.lastChrom = this.curChr;
	}

	
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TabularOutput#headerStr()
	 */
	public String headerStr() throws Exception {
		return null;
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TabularOutput#toCsvStr()
	 */
	public String toCsvStr() throws Exception {
		return null;
	}

	
}
