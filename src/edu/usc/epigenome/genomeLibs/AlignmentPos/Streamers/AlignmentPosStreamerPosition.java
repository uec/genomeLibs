/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Symbol;

import BisulfiteCytosines.CpgPair;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.Counters.NmerCounter;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicPositionScored;

/**
 * @author benb
 *
 */
public class AlignmentPosStreamerPosition {

	/**
	 * @param preWindSize
	 * @param postWindSize
	 */
	public AlignmentPosStreamerPosition(int preWindSize, int postWindSize) {
		super();
		this.preWindSize = preWindSize;
		this.postWindSize = postWindSize;
	}


	public int preWindSize = 0;
	public int postWindSize = 0;
	
	public AlignmentPos[] priorAps = null;
	public AlignmentPos currentAp = null;
	public GenomicPositionScored currentScoredPosition = null; // By default , this is the AP itself.  But it can be overridden.
	public AlignmentPos[] nextAps = null;
	
	public NmerCounter preNmerCounts = null;
	public NmerCounter nextNmerCounts = null;
	
	public CpgPair cpg = null; // Does not need to be supplied.
	
	public static final Symbol[] CPG_SYMBOLS = { DNATools.c(), DNATools.g() };
	
	public double preNmerCpgDensity()
	{
		double out = 0.0;
		if (preNmerCounts!=null)
		{
			int nCpg = preNmerCounts.getCount(CPG_SYMBOLS);
			out = (double)nCpg / (double)priorAps.length;
		}
		return out;
	}
	
	public double nextNmerCpgDensity()
	{
		double out = 0.0;
		if (nextNmerCounts!=null)
		{
			int nCpg = nextNmerCounts.getCount(CPG_SYMBOLS);
			out = (double)nCpg / (double)nextAps.length;
		}
		return out;
	}	
	
	public int getWindSize()
	{
		return preWindSize + postWindSize + 1;
	}
	
	public double getAvgScore()
	{
			
		double total = 0.0;
		for (AlignmentPos ap : priorAps)
		{
			total += ap.getSummaryScore();
		}
		total += currentAp.getSummaryScore();
		for (AlignmentPos ap : nextAps)
		{
			total += ap.getSummaryScore();
		}
		
		//System.err.println("total = " + total + "\tWind size = " +this.getWindSize() );
		return total / (double)this.getWindSize();
	}
	
	public double getAvgDepth()
	{
			
		double total = 0.0;
		for (AlignmentPos ap : priorAps)
		{
			total += ap.getTotalDepth();
		}
		total += currentAp.getTotalDepth();
		for (AlignmentPos ap : nextAps)
		{
			total += ap.getTotalDepth();
		}
		
		return total / (double)this.getWindSize();
	}
	
}
