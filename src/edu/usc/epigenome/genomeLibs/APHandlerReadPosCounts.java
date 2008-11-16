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
public class APHandlerReadPosCounts extends ReadPosCounts implements AlignmentPosStreamHandler {

	/**
	 * 
	 */
	public APHandlerReadPosCounts() {
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	public boolean streamElement(AlignmentPos[] priorAps,
			AlignmentPos currentAp, AlignmentPos[] nextAps) 
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
