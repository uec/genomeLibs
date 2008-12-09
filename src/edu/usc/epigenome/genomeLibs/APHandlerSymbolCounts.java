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
public class APHandlerSymbolCounts extends ReadPosCounter implements AlignmentPosStreamHandler {

	/**
	 * 
	 */
	private static final long serialVersionUID = -3266222243481309526L;



	/**
	 * 
	 */
	public APHandlerSymbolCounts() {
	}

	

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	public boolean streamElement(AlignmentPos[] pre, AlignmentPos currentAp, AlignmentPos[] post) 
	{
		boolean passes = true;
		Iterator<ReadPos> it = currentAp.getReadPositions().iterator();
		while (it.hasNext())
		{
			ReadPos rp = it.next();
			this.increment(rp);
		}
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
