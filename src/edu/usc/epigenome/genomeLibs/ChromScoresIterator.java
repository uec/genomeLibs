package edu.usc.epigenome.genomeLibs;

import java.util.Iterator;


// Returns non-overlapping sections of the genome, usually
// entire chromosomes.
abstract public class ChromScoresIterator implements Iterator<ChromScoresFast> {

	
	/***********
	 * Interface functions
	 */
	
	abstract public boolean hasNext();
	abstract public ChromScoresFast next();

	
	public void remove()
	throws UnsupportedOperationException 
	{
		throw new UnsupportedOperationException();
	}
	
	
	
}
