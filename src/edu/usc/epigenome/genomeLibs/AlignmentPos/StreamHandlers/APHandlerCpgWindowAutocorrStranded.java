/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

/**
 * @author benb
 * 
 *
 */
 public class APHandlerCpgWindowAutocorrStranded extends APHandlerCpgWindowStats {

	 public int[] totalsSame = null;
	 public double[] totalDiffsSame = null;
	 public int[] totalsOpposite = null;
	 public double[] totalDiffsOpposite = null;

	 public int uniqueCpgs = 0;	
	 public int uniqueCpgsSame = 0;	
	 public int uniqueCpgsOpposite = 0;	
	 

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	/**
	 * @param inWindSize
	 */
	public APHandlerCpgWindowAutocorrStranded(int inWindSize) {
		super(inWindSize);
	}

	public void init() {
		super.init();
		
		totalsSame = new int[this.windSize];
		totalDiffsSame = new double[this.windSize];
		totalsOpposite= new int[this.windSize];
		totalDiffsOpposite = new double[this.windSize];

	}

	public void finish() {
		super.finish();

		ListUtils.setDelim(",");
		System.out.println(ListUtils.excelLine(totalsSame));
		System.out.println(ListUtils.excelLine(totalDiffsSame));
		System.out.println(ListUtils.excelLine(totalsOpposite));
		System.out.println(ListUtils.excelLine(totalDiffsOpposite));

		System.err.println("Unique CpGs contributing: " + uniqueCpgs);
		System.err.println("Unique CpGs contributing (same strand): " + uniqueCpgsSame);
		System.err.println("Unique CpGs contributing (opposite strand): " + uniqueCpgsOpposite);

	
	}


	public boolean streamWindow(AlignmentPosStreamerPosition streamPos, CpgPair pair, Queue<CpgPair> cpgWind)
	{
		int thisPos = pair.getPos();
		double thisMethFw = pair.getMethylatedFrac(true);
		double thisMethRev = pair.getMethylatedFrac(false);

		boolean countIt = false; // For counting num Cpgs used
		boolean countItSame = false; // For counting num Cpgs used
		boolean countItOpposite = false; // For counting num Cpgs used
		
//		System.err.println("Starting " + thisPos + " --------------");
		
		// Add our self count
		if (!Double.isNaN(thisMethFw) && !Double.isNaN(thisMethRev))
		{
			totalsOpposite[0]++;
			totalDiffsOpposite[0] += Math.abs(thisMethFw-thisMethRev);
			countIt = true;
			countItOpposite = true;
//			System.err.println(thisPos + " -> " + thisPos + "\tself\t" + Math.abs(thisMethFw-thisMethRev));
		}
		
		
		for( CpgPair otherPair : cpgWind)
		{
			int otherPos = otherPair.getPos();
			int dist = thisPos - otherPos;

			double otherMethFw = otherPair.getMethylatedFrac(true);
			double otherMethRev = otherPair.getMethylatedFrac(false);
			
			
			// FW vs. FW
			if (!Double.isNaN(thisMethFw) && !Double.isNaN(otherMethFw))
			{
				totalsSame[dist]++;
				totalDiffsSame[dist] += Math.abs(thisMethFw-otherMethFw);
				countIt = true;
				countItSame = true;
//				System.err.println(otherPos + " -> " + thisPos + "\tFW/FW\t" + Math.abs(thisMethFw-otherMethFw));
			}
			
			
			// REV vs. REV
			if (!Double.isNaN(thisMethRev) && !Double.isNaN(otherMethRev))
			{
				totalsSame[dist]++;
				totalDiffsSame[dist] += Math.abs(thisMethRev-otherMethRev);
				countIt = true;
				countItSame = true;
//				System.err.println(otherPos + " -> " + thisPos + "\tREV/REV\t" + Math.abs(thisMethRev-otherMethRev));
			}
			
			
			// FW vs. REV
			if (!Double.isNaN(thisMethFw) && !Double.isNaN(otherMethRev))
			{
				totalsOpposite[dist]++;
				totalDiffsOpposite[dist] += Math.abs(thisMethFw-otherMethRev);
				countIt = true;
				countItOpposite = true;
//				System.err.println(otherPos + " -> " + thisPos + "\tFW/REV\t" + Math.abs(thisMethFw-otherMethRev));
			}			
			
			// REV vs. FW
			if (!Double.isNaN(thisMethRev) && !Double.isNaN(otherMethFw))
			{
				totalsOpposite[dist]++;
				totalDiffsOpposite[dist] += Math.abs(thisMethRev-otherMethFw);
				countIt = true;
				countItOpposite = true;
//				System.err.println(otherPos + " -> " + thisPos + "\tREV/FW\t" + Math.abs(thisMethRev-otherMethFw));
			}			

		}

		if (countIt) uniqueCpgs++;
		if (countItSame) uniqueCpgsSame++;
		if (countItOpposite) uniqueCpgsOpposite++;
		
		return countIt;
	}


}
