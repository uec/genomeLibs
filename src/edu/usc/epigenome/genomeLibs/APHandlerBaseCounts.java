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
 */
public class APHandlerBaseCounts extends AlignmentPosStreamHandler {

	public int MAX_CYCLES = 100;
	
	protected HashMap<ReadPos,Integer> cycleCounts;
	
	/**
	 * 
	 */
	public APHandlerBaseCounts() {
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#finish()
	 */
	@Override
	public void finish() {
		// Nothing to do
		System.err.println("Finishing APHandlerBaseCounts");
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#init()
	 */
	@Override
	public void init() {
		// Initialize maps
		System.err.println("Initializing APHandlerBaseCounts");
		cycleCounts = new HashMap<ReadPos,Integer>();
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	@Override
	public boolean streamElement(LinkedList<AlignmentPos> priorAps,
			AlignmentPos currentAp, LinkedList<AlignmentPos> nextAps) 
	{
		boolean passes = true;
		
		AlignmentPosSnps currentApSnps = null;
		try
		{
			currentApSnps = (AlignmentPosSnps)currentAp;
		}
		catch (ClassCastException e)
		{
			System.err.println("APHandlerBaseCounts called with non-SNP AlignmentPos objects");
			e.printStackTrace();
			System.exit(0);
		}
		
		Iterator<ReadPos> rpIt = currentApSnps.getReadPositions().iterator();
		while (rpIt.hasNext())
		{
			// ReadPos has equals re-implemented so that identical ones are equal
			ReadPos rp = rpIt.next();
			Integer count = cycleCounts.get(rp);
			if (count == null) count = new Integer(0);
			count++;
			cycleCounts.put(rp, count);
		}
		
		return passes;
	}

	
	
	
	
	/******  OUTPUT ********/
	
	
	public String excelOutput()
	throws IllegalSymbolException
	{
		String out = "";
		
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
