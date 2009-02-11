package edu.usc.epigenome.genomeLibs;

import java.util.Collection;

import org.biojava.bio.seq.DNATools;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosOptions;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosSnpsBisulfiteConverted;
/**
 * 
 */

/**
 * @author benb
 *
 * Holds a pair of CpG AlignmentPos objects, and knows what to do with them.  They should be 
 * flipped so that both are relative to the cytosine.
 */
public class CpgPair {

	protected AlignmentPosSnpsBisulfiteConverted cpgF = null;
	protected AlignmentPosSnpsBisulfiteConverted cpgR = null;

	/**
	 * @param cpgA
	 * @param cpgB
	 */
	public CpgPair(AlignmentPosSnpsBisulfiteConverted cpgA,
			AlignmentPosSnpsBisulfiteConverted cpgB) {
		super();
		this.cpgF = cpgA;
		this.cpgR = cpgB;
		
		if (cpgA.getRef() != DNATools.c()) System.err.println("CpgPair: Why is CpgF not a cytosine?");
		if (cpgB.getRef() != DNATools.c()) System.err.println("CpgPair: Why is CpgR not a cytosine?");
	}
	
	public String getChr()
	{
		return cpgF.getChr();
	}
	
	public int getPos()
	{
		return cpgF.getPos();
	}
	
	public int getDepth(boolean fw)
	{
		AlignmentPosSnpsBisulfiteConverted c = (fw) ? cpgF : cpgR;
		return c.getDepth(true);
	}

	public int getDepthWithIdentical(boolean fw)
	{
		AlignmentPosSnpsBisulfiteConverted c = (fw) ? cpgF : cpgR;
		return c.getDepthWithIdentical(true);
	}

	public int getDepthNoIdentical(boolean fw)
	{
		AlignmentPosSnpsBisulfiteConverted c = (fw) ? cpgF : cpgR;
		return c.getDepthNoIdentical(true);
	}

	public double getMethylatedFrac(boolean fw)
	{
		AlignmentPosSnpsBisulfiteConverted c = (fw) ? cpgF : cpgR;
		return c.getMethylatedFrac();
	}
	
	public double getMethylatedFrac(boolean fw, int inMaxIdentical)
	{
		AlignmentPosSnpsBisulfiteConverted c = (fw) ? cpgF : cpgR;
		return c.getMethylatedFrac(inMaxIdentical);
	}

	public double getStrandMismatch()
	{
		double hemi = Math.abs(getMethylatedFrac(true) - getMethylatedFrac(false));
		return hemi;
	}
	
	public String getMethylatedFracStr(boolean fw)
	{
		AlignmentPosSnpsBisulfiteConverted c = (fw) ? cpgF : cpgR;
		return c.getMethylatedFracString();
	}

	public String getMethylatedFracStr(boolean fw, int inMaxIdentical)
	{
		AlignmentPosSnpsBisulfiteConverted c = (fw) ? cpgF : cpgR;
		return c.getMethylatedFracString(inMaxIdentical);
	}

	
	public String gffLine()
	{
		String out = "";
		
		out += cpgF.getChr();
		
		out += "\tCpG";
		out += "\texon";
		out += "\t" + cpgF.getPos();
		out += "\t" + cpgR.getPos();
		out += "\t" + this.getMethylatedFracStr(true);
		out += "\t" + cpgF.getStrandSymbol();
		out += "\t.";
		
		// Group sec
		out += "\t";
		
		out += "fw_depth " + this.getDepth(true) + "; ";
		out += "rev_depth " + this.getDepth(false) + "; ";

		out += "fw_conv " + this.getMethylatedFracStr(true) + "; ";
		out += "rev_conv " + this.getMethylatedFracStr(false) + "; ";

		out += "strandMismatch " + this.getStrandMismatch() + "; ";
	
//		out += cpgA.gffLine();
//		out += cpgB.gffLine();
		
		return out;
	}
	
	public String csvLine()
	{
		String out = "";
		
		out += cpgF.getChr();
		
		out += "," + cpgF.getPos();
		out += "," + cpgR.getPos();
		out += "," + cpgF.getStrandSymbol();
		out += "," + this.csvStats();
		
		return out;
	}
	
	public String csvStats()
	{
		double[] myStats = this.stats();
		return csvStats(myStats);
	}
	
	public static String csvStats(Collection<CpgPair> cpgs)
	{
		double[] myStats = stats(cpgs);
		return csvStats(myStats);
	}
	
	protected static String csvStats(double[] inStats)
	{
		ListUtils.setDelim(",");
		String out = ListUtils.excelLine(inStats);
		return out;
	}
		
	public double[] stats()
	{
		double[] out = new double[9];

		out[0] = (double)this.getDepthWithIdentical(true);
		out[1] = (double)this.getDepthNoIdentical(true);
		out[2] = (double)this.getMethylatedFrac(true, 0);
		out[3] = (double)this.getMethylatedFrac(true, 1);
		out[4] = (double)this.getDepthWithIdentical(false);
		out[5] = (double)this.getDepthNoIdentical(false);
		out[6] = (double)this.getMethylatedFrac(false, 0);
		out[7] = (double)this.getMethylatedFrac(false, 1);
		out[8] = (double)this.getStrandMismatch();
		
		return out;
	}
	
	public static double[] emptyStats()
	{
		double[] out = new double[9];

		out[0] = Double.NaN;
		out[1] = Double.NaN;
		out[2] = Double.NaN;
		out[3] = Double.NaN;
		out[4] = Double.NaN;
		out[5] = Double.NaN;
		out[6] = Double.NaN;
		out[7] = Double.NaN;
		out[8] = Double.NaN;
		
		return out;
	}

	public static double[] stats(Collection<CpgPair> cpgs)
	{
		double[] totals = null;
		double[] ns = null; // So we can have NaNs
		
		for ( CpgPair cpg : cpgs)
		{
			double[] cpgStats = cpg.stats();

			// Create the output array if we haven't already
			if (totals==null)
			{
				totals = new double[cpgStats.length];
				ns = new double[cpgStats.length];
			}
			
			// Total only those that are not NAN
			for (int i = 0; i < cpgStats.length; i++)
			{
				if (!Double.isInfinite(cpgStats[i]) && !Double.isNaN(cpgStats[i]))
				{
					totals[i] += cpgStats[i];
					ns[i]++;
				}
			}
		}
		
		double[] out = null;
		if (totals==null)
		{
			out = emptyStats();
		}
		else
		{
			out = new double[totals.length];
			for (int i = 0; i < totals.length; i++)
			{
				out[i] = totals[i] / ns[i]; // NaNs are ok.
			}
		}
		
		return out;
	}
	
	
}
