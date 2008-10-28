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
	
	protected HashMap<Symbol,Integer>[] cycleCounts = ((HashMap<Symbol,Integer>[])new HashMap[MAX_CYCLES]);
	
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
		for (int i = 0; i < MAX_CYCLES; i++)
		{
			cycleCounts[i] = new HashMap<Symbol,Integer>();
		}
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
			ReadPos rp = rpIt.next();
			Symbol sym = rp.getSymReaddir();
			int cycle = rp.getReadPos();
			Integer count = cycleCounts[cycle].get(sym);
			if (count == null)
			{
				count = new Integer(0);
			}
			count++; // Does this actually modify the stored object?
		}
		
		return passes;
	}

	
	
	
	
	/******  OUTPUT ********/
	
	public int[][] countMatrix()
	{
		int[][] out = new int[MAX_CYCLES][DNATools.getDNA().size()];
		for (int i = 0; i < MAX_CYCLES; i++)
		{
			HashMap<Symbol,Integer> map = cycleCounts[i];

			Iterator<Symbol> it = ((Iterator<Symbol>)DNATools.getDNA().iterator());
			int symInd = 0;
			while(it.hasNext())
			{
				Integer count = map.get(it.next());
				if (count == null) count = new Integer(0);
				out[i][symInd] = count.intValue();
				symInd++;
			}
		}
		
		// Trim to max cycle and max symbol size
		
		return out;
	}
	
	public String excelOutput()
	throws IllegalSymbolException
	{
		String out = "";
		
		// Header
		List<Object> l = new Vector<Object>();
		Iterator<Symbol> it = (Iterator<Symbol>)DNATools.getDNA().iterator();
		while (it.hasNext())
		{
			l.add(new Character(DNATools.dnaToken(it.next())));
		}
		out += ListUtils.tabbedLine(l) + "\n";
		
		// And the counts
		int[][] mat = this.countMatrix();
		out += MatUtils.matString(mat);
		
		return out;
	}

}
