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
public class APHandlerSymbolCounts extends SymbolCounts implements AlignmentPosStreamHandler {

	/**
	 * 
	 */
	public APHandlerSymbolCounts() {
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
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

	
	
	

}
