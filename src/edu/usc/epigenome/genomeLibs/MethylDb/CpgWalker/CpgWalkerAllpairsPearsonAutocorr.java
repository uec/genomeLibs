package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import com.mallardsoft.tuple.Pair;
import com.mallardsoft.tuple.Quadruple;
import com.mallardsoft.tuple.Tuple;

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

public class CpgWalkerAllpairsPearsonAutocorr extends CpgWalkerAllpairs {

	protected double methMean = Double.NaN;
	protected double methSD = Double.NaN;

	protected double[] counts;
	protected double[] sums;
	
	
	public CpgWalkerAllpairsPearsonAutocorr(CpgWalkerParams inWalkParams,
			boolean inSamestrandOnly, double inMethMean, double inMethStandardDev) {
		super(inWalkParams, inSamestrandOnly);
		
		methMean = inMethMean;
		methSD = inMethStandardDev;
		
		init();
	}
	
	
	protected void init()
	{
		useSummarizers = false;
		
		// Initalize counters
		int windSize = this.walkParams.maxWindSize;
		counts = new double[windSize-1]; 
		sums = new double[windSize-1];
	}

	
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#reset()
	 */
	@Override
	public void reset() {
		super.reset();
		init();
	}


	@Override
	protected void recordPair(Cpg first, Cpg second)
	{
		// Get the distance
		
		int dist = second.chromPos - first.chromPos - 1; // ([0] will be dist of 1)
		double m1 = first.fracMeth(true);
		double m2 = second.fracMeth(true);
		
		if (Double.isNaN(methMean) || Double.isNaN(methSD))
		{
			System.err.println("mean or SD must be set before running CpgWalkerAllpairsPearsonAutocorr::recordPair()");
			System.exit(1);
		}
		
		// (only count it if both have a valid range)
		if (!Double.isNaN(m1) && !Double.isNaN(m2)) 
		{
			this.counts[dist]++;
			
			double p1 = (m1 - this.methMean) / this.methSD;
			double p2 = (m2 - this.methMean) / this.methSD;
			this.sums[dist] += (p1 * p2);
		}
		
	}
	
//	@Override
//	public String headerStr()
//	{
//		StringBuffer sb = new StringBuffer((int)1E5);
//		
//		sb.append("Range");
//		for (int i = 0; i < this.walkParams.maxWindSize; i++)
//		{
//			sb.append(",d");
//			sb.append(i);
//		}
//		
//		return sb.toString();
//	}
	
	

	
	@Override
	public String toCsvStr()
	throws Exception
	{
		StringBuffer sb = new StringBuffer((int)1E6);
		ListUtils.setDelim(",");

		double[] countsMinusOne = MatUtils.sumVect(this.counts, -1.0);
		// Values of -1.0 will now give us a pearson of 0 if we don't set them back to 0
		for (int i = 0; i < counts.length; i++)
		{
			if (countsMinusOne[i]==-1.0) countsMinusOne[i] = 0.0; 
		}
		double[] pearsons = MatUtils.divVects(this.sums, countsMinusOne); 
		
		sb.append(ListUtils.excelLine(pearsons));
		sb.append(ListUtils.excelLine(countsMinusOne));
		
		return sb.toString();
	}


	

}
