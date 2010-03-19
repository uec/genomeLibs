package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import edu.usc.epigenome.genomeLibs.TabularOutput;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

public abstract class CpgWalkerMultisample implements TabularOutput,ActionListener {

	public CpgWalkerParams walkParams = null;
	
	protected CpgWalker[] walkers = null;
	
	public CpgWalkerMultisample(CpgWalkerParams inWalkParams, int numSamples) {
		this.walkParams=inWalkParams;
		this.walkers = new CpgWalker[numSamples];
		for (int i = 0; i < numSamples; i++)
		{
			walkers[i] = new CpgWalker(walkParams);
		}
		
		// And add ourselves as a listener.  Add it to the
		// last one, so we know that all of them have been
		// processed.
		walkers[numSamples-1].addWindowListener(this);
	}

	public int numSamples()
	{
		return walkers.length;
	}
	
	public void newChrom()
	{
		reset();
	}
	
	public String getCurChr()
	{
		return walkers[0].getCurChr();
	}

	public void setCurChr(String inChr)
	{
		for (CpgWalker w : walkers)
		{
			w.setCurChr(inChr);
		}
		reset();
	}
	
	public void finishChr()
	{
		for (CpgWalker w : walkers)
		{
			w.finishChr();
		}
	}

	public void reset()
	{
		for (CpgWalker w : walkers)
		{
			w.reset();
		}
	}
	
	public void resetSummarizers()
	{
		for (CpgWalker w : walkers)
		{
			w.resetSummarizers();
		}
	}	
	
	public void streamCpg(Cpg[] cpgs)
	{
		for (int i = 0; i < numSamples(); i++)
		{
			CpgWalker w = walkers[i];
			w.streamCpg(cpgs[i]);
		}
	}
	
	public int windLen()
	{
		return walkers[0].windLen();
	}
	
	public int windStart()
	{
		return walkers[0].windStart();
	}

	public int windEnd()
	{
		return walkers[0].windLen();
	}
	
	public String lastProcessedWindStr(boolean longVers)
	{
		return this.walkers[0].lastProcessedWindStr(longVers);
	}

	public String windStr(boolean longVers)
	{
		return this.walkers[0].windStr(longVers);
	}

	public String windStr()
	{
		return this.windStr(false);
	}
	
	public void actionPerformed(ActionEvent e) {
		// we should get an event when our last CpgWindow
		// gets a good window.
		if (e.getID() == CpgWalker.PROCESS_WINDOW_EVENT)
		{
			this.processWindow();
		}
		else
		{
			System.err.println("CpgWalkerMultisample received illegal event " + e.getID());
		}
		
	}
	
	// **** OVERRIDE THIS ******
	abstract public void processWindow();

	public String headerStr() throws Exception {
		// TODO Auto-generated method stub
		return null;
	}

	public String toCsvStr() throws Exception {
		// TODO Auto-generated method stub
		return null;
	}

}
