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
public class APHandlerSNPCounts extends StringCounter implements AlignmentPosStreamHandler {

	/**
	 * Constructor
	 */
	public APHandlerSNPCounts() {
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
			String key = AlignmentPos.getRefTokens(priorAps) +
			"/" + 
			currentAp.getRefToken(true) + 
			"/" + 
			AlignmentPos.getRefTokens(nextAps) + 
			" -> " +
			rp.getSymToken() + rp.getStrandChar();
			
//			System.err.println("\t" + key);
			
			this.increment(key);
		}
		
		return passes;
	}

	
	
	

}
