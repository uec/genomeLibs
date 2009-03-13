/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import BisulfiteCytosines.CpgPair;

import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

/**
 * @author benb
 * 
 *
 */
 public class APHandlerCpgWindowAutocorr extends APHandlerCpgWindowStats {

	 public int[] totals = null;
	 public double[] totalDiffs = null;
		

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	/**
	 * @param inWindSize
	 */
	public APHandlerCpgWindowAutocorr(int inWindSize) {
		super(inWindSize);
	}

	public void init() {
		super.init();
		
		totals = new int[this.windSize];
		totalDiffs = new double[this.windSize];
		
//		for (int i = 0; i < this.windSize; i++)
//		{
//			totals[i] = -1;
//			totalDiffs[i] = -1.0;
//		}
	}

	public void finish() {
		super.finish();

		ListUtils.setDelim(",");
		System.out.println(ListUtils.excelLine(totals));
		System.out.println(ListUtils.excelLine(totalDiffs));
	}


	public boolean streamWindow(AlignmentPosStreamerPosition streamPos, CpgPair pair, Queue<CpgPair> cpgWind)
	{
		int thisPos = pair.getPos();
		double thisMeth = pair.getMethylatedFrac();
		
		for( CpgPair otherPair : cpgWind)
		{
			int dist = thisPos - otherPair.getPos();
			double diff = Math.abs(thisMeth - otherPair.getMethylatedFrac());

			totals[dist-1]++;
			totalDiffs[dist-1] += diff;
		}

		
		return true;
	}


}
