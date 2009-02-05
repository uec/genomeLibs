/**
 * 
 */
package edu.usc.epigenome.genomeLibs.ReadPos.StreamHandlers;

import org.biojava.bio.symbol.Symbol;

import edu.usc.epigenome.genomeLibs.Counters.SymbolCounterStratifiedByCycle;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;

/**
 * @author benb
 * 
 *
 */
public class RPHandlerSymbolCountsStratifyByCycle extends SymbolCounterStratifiedByCycle implements ReadPosStreamHandler {

	private static final long serialVersionUID = 5505840791821073036L;


	/**
	 * Constructor
	 */
	public RPHandlerSymbolCountsStratifyByCycle() {
		super();
	}

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	public void init() {
	}

	public void finish() {
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	public boolean streamElement(ReadPos currentRp) 
	{
		boolean passes = true;
		Symbol sym = currentRp.getSym();
		int cyc = currentRp.getCycle();

		if ((cyc>=getLowCycleStart())&&(cyc<=getLowCycleEnd()))
		{
			this.increment(sym, getLowStratString());
		}
		else if ((cyc>=getMediumCycleStart())&&(cyc<=getMediumCycleEnd()))
		{
			this.increment(sym, getMediumStratString());
		}
		else if ((cyc>=getHighCycleStart())&&(cyc<=getHighCycleEnd()))
		{
			this.increment(sym, getHighStratString());
		}

		return passes;
	}

	

}
