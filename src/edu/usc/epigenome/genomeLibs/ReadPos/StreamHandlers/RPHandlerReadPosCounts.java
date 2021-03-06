/**
 * 
 */
package edu.usc.epigenome.genomeLibs.ReadPos.StreamHandlers;

import edu.usc.epigenome.genomeLibs.Counters.ReadPosCounter;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;

/**
 * @author benb
 * 
 *
 */
public class RPHandlerReadPosCounts extends ReadPosCounter implements ReadPosStreamHandler {

	private static final long serialVersionUID = -8766957148385079068L;



	/**
	 * Constructor
	 */
	public RPHandlerReadPosCounts() {
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
		this.increment(currentRp);
		return passes;
	}



	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TreeMapCounter#excelOutput()
	 */
	@Override
	public String excelOutput() {
		// TODO Auto-generated method stub
		return super.excelOutput();
	}

	
	
	

}
