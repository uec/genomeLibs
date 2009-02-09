package edu.usc.epigenome.genomeLibs;

import org.biojava.bio.seq.DNATools;

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
	
	public int getDepth(boolean fw)
	{
		AlignmentPosSnpsBisulfiteConverted c = (fw) ? cpgF : cpgR;
		return c.getDepth(true);
	}

	public double getMethylatedFrac(boolean fw)
	{
		AlignmentPosSnpsBisulfiteConverted c = (fw) ? cpgF : cpgR;
		return c.getMethylatedFrac();
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
		out += "," + this.getMethylatedFracStr(true);
		out += "," + cpgF.getStrandSymbol();
		out += "," + this.getDepth(true);
		out += "," + this.getDepth(false);
		out += "," + this.getMethylatedFracStr(true);
		out += "," + this.getMethylatedFracStr(false);
		out += "," + this.getStrandMismatch();
		
		return out;
	}
		
	
}
