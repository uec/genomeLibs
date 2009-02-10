/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Symbol;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.Counters.NmerCounter;

/**
 * @author benb
 *
 */
public class AlignmentPosStreamerPosition {

	public AlignmentPos[] priorAps = null;
	public AlignmentPos currentAp = null;
	public AlignmentPos[] nextAps = null;
	
	public NmerCounter preNmerCounts = null;
	public NmerCounter nextNmerCounts = null;
	
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
	
}
