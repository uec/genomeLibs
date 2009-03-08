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
		
//		if (cpgA.getRef() != DNATools.c()) System.err.println("CpgPair: Why is CpgF not a cytosine?");
//		if (cpgB.getRef() != DNATools.c()) System.err.println("CpgPair: Why is CpgR not a cytosine?");
	}
	
	public String getChr()
	{
		return cpgF.getChr();
	}
	
	public int getPos()
	{
		return cpgF.getPos();
	}
	
	
	
	// *******************************************
	//
	//   DEPTH
	// 
	// *******************************************
	
	
	
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

	// *******************************************
	//
	//   SNP
	// 
	// *******************************************
	
	public double getSnpProbCombined()
	{
		double fw = getSnpProb(true);
		double rev = getSnpProb(false);
		
		return Math.max(fw, rev);
	}

	public double getSnpProb(boolean fw)
	{
		AlignmentPosSnpsBisulfiteConverted c = (fw) ? cpgF : cpgR;
		return c.getSnpProb();
	}
	
	public double getSnpProbCombined(int inMaxIdentical)
	{
		double fw = getSnpProb(true, inMaxIdentical);
		double rev = getSnpProb(false, inMaxIdentical);
		
		return Math.max(fw, rev);
	}

	public double getSnpProb(boolean fw, int inMaxIdentical)
	{
		AlignmentPosSnpsBisulfiteConverted c = (fw) ? cpgF : cpgR;
		return c.getSnpProb(inMaxIdentical);
	}
	
	// *******************************************
	//
	//   Methylation
	// 
	// *******************************************

	public double getMethylatedFrac()
	{
		double fwDepth = (double)(cpgF.getDepth(true));
		double revDepth = (double)(cpgR.getDepth(true));

		double fwFrac = (fwDepth>0) ? (fwDepth/(fwDepth+revDepth)) : 0;
		double revFrac = 1-fwFrac;
		
		double fwMeth = cpgF.getMethylatedFrac();
		double revMeth = cpgR.getMethylatedFrac();
		
		// If it's a SNP, "depth" can be 0, but the number of C's and T's can be 0.
		// In this case, fwMeth will be Nan
		double total = 0.0;
		total += ((fwDepth > 0) && !Double.isNaN(fwMeth)) ? (fwFrac*fwMeth) : 0.0;
		total += ((revDepth > 0) && !Double.isNaN(revMeth)) ? (revFrac*revMeth) : 0.0;
		
		return total;
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

	// *******************************************
	//
	//   GFF
	// 
	// *******************************************

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
	
	// *******************************************
	//
	//   CSV STATS
	// 
	// *******************************************
	
	
	
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
	
	public static String csvStatsHeaders()
	{
		ListUtils.setDelim(",");
		String labels[] = {
				// FW
				"fwDepth", 
				"fwDepthCollapsed", 
				"fwMeth", 
				"fwMethCollapsed", 
				"fwSnpProb", 
				"fwSnpProbCollapsed", 
				//REV
				 "revDepth", 
				 "revDepthCollapsed", 
				 "revMeth", 
				 "revMethCollapsed",
				 "revSnpProb", 
				 "revSnpProbCollapsed", 
				 //OTHER
				 "methOverall",
				 "strandMismatch" 
				 }; 
		String out = ListUtils.excelLine(labels);
		return out;
	}	
	
	public double[] stats()
	{
		double[] out = new double[14];

		int ind = 0;
		
		// FW
		out[ind++] = (double)this.getDepthWithIdentical(true);
		out[ind++] = (double)this.getDepthNoIdentical(true);
		out[ind++] = (double)this.getMethylatedFrac(true, 0);
		out[ind++] = (double)this.getMethylatedFrac(true, 1);
		out[ind++] = (double)this.getSnpProb(true, 0);
		out[ind++] = (double)this.getSnpProb(true, 1);
		
		// REV
		out[ind++] = (double)this.getDepthWithIdentical(false);
		out[ind++] = (double)this.getDepthNoIdentical(false);
		out[ind++] = (double)this.getMethylatedFrac(false, 0);
		out[ind++] = (double)this.getMethylatedFrac(false, 1);
		out[ind++] = (double)this.getSnpProb(false, 0);
		out[ind++] = (double)this.getSnpProb(false, 1);
		
		// OTHER
		out[ind++] = (double)this.getMethylatedFrac();
		out[ind++] = (double)this.getStrandMismatch();
		
		return out;
	}
	
	public static double[] emptyStats()
	{
		int nStats = 14;
		
		double[] out = new double[14];
		for (int i = 0; i < nStats; i++)
		{
			out[i] = Double.NaN;
		}
		
		return out;
	}

	public static double[] stats(Collection<CpgPair> cpgs)
	{
		double[] totals = null;
		double[] ns = null; // So we can have NaNs
		
		for ( CpgPair cpg : cpgs)
		{
			double[] cpgStats = cpg.stats();
			//System.err.println("Window has " + cpgs.size() + " cpgs");

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
