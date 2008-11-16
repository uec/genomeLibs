/**
 * 
 */
package edu.usc.epigenome.genomeLibs;

import java.util.*;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

/**
 * @author benb
 * 
 *
 */
public class StringCounts  {

	protected TreeMap<String,Integer> counts;
	
	/**
	 * 
	 */
	public StringCounts() {
	}
	
	public StringCounts(StringCounts inSc)
	{
		counts = (TreeMap<String,Integer>)inSc.counts.clone();
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#finish()
	 */
	
	public void finish() {
		// Nothing to do
		// System.err.println("Finishing SymbolCounts");
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#init()
	 */
	public void init() {
		// Initialize maps
		//System.err.println("Initializing SymbolCounts");
		counts = new TreeMap<String,Integer>();
	}

	
	
	/******  OUTPUT ********/
	
	
	public String excelOutput()
	throws IllegalSymbolException
	{
		String out = "";
		
		// System.err.println("rpIt has " + this.cycleCounts.keySet().size() + " elements");
		Iterator<String> strIt = this.counts.keySet().iterator();
		while (strIt.hasNext())
		{
			String symbol = strIt.next();
			out += symbol + "," + this.counts.get(symbol);
			out += "\n";
		}
		
		return out;
	}

}
