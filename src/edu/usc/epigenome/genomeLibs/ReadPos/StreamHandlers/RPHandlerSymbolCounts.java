/**
 * 
 */
package edu.usc.epigenome.genomeLibs.ReadPos.StreamHandlers;

import edu.usc.epigenome.genomeLibs.Counters.SymbolCounter;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;

/**
 * @author benb
 * 
 *
 */
public class RPHandlerSymbolCounts extends SymbolCounter implements ReadPosStreamHandler {

	private static final long serialVersionUID = 2793815144706813135L;



	/**
	 * Constructor
	 */
	public RPHandlerSymbolCounts() {
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
		this.increment(currentRp.getSym());
		return passes;
	}



	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TreeMapCounter#excelOutput()
	 */
	@Override
	public String excelOutput() {
		return super.excelOutput();
	}

	
	
	

}
