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
public class SymbolCounts  {

	protected HashMap<ReadPos,Integer> cycleCounts;
	
	/**
	 * 
	 */
	public SymbolCounts() {
	}
	
	public SymbolCounts(SymbolCounts inSc)
	{
		cycleCounts = (HashMap<ReadPos,Integer>)inSc.cycleCounts.clone();
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#finish()
	 */
	
	public void finish() {
		// Nothing to do
		System.err.println("Finishing SymbolCounts");
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#init()
	 */
	public void init() {
		// Initialize maps
		System.err.println("Initializing SymbolCounts");
		cycleCounts = new HashMap<ReadPos,Integer>();
	}

	
	
	/******  OUTPUT ********/
	
	
	public String excelOutput()
	throws IllegalSymbolException
	{
		String out = "";
		
		// System.err.println("rpIt has " + this.cycleCounts.keySet().size() + " elements");
		Iterator<ReadPos> rpIt = this.cycleCounts.keySet().iterator();
		while (rpIt.hasNext())
		{
			ReadPos rp = rpIt.next();
			out += rp.commaSeparatedLine();
			out += "," + this.cycleCounts.get(rp);
			out += "\n";
		}
		
		return out;
	}

}
