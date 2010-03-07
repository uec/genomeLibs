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

	protected double methMeans[] = null;
	protected double methSDs[] = null;

	protected double[] counts;
	protected double[] pearsonSums;
	
	// These are for calculating mean and sd
	protected double[] meanCounts;
	protected double[] methTotals;
	protected double[] methSquareTotals;
	
	public CpgWalkerAllpairsPearsonAutocorr(CpgWalkerParams inWalkParams,
			boolean inSamestrandOnly, double[] inMethMeans, double inMethStandardDevs[]) {
		super(inWalkParams, inSamestrandOnly);
		
		methMeans = inMethMeans;
		methSDs = inMethStandardDevs;
		
		init();
	}
	
	
	protected void init()
	{
		useSummarizers = true;
		
		// Initalize counters
		int windSize = this.walkParams.maxScanningWindSize;
		counts = new double[windSize-1]; 
		meanCounts = new double[windSize-1]; 
		pearsonSums = new double[windSize-1];
		methTotals = new double[windSize-1];
		methSquareTotals = new double[windSize-1];
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
		
		if ((methMeans==null) || (methSDs==null))
		{
			System.err.println("mean or SD must be set before running CpgWalkerAllpairsPearsonAutocorr::recordPair()");
			System.exit(1);
		}
		
		// (only count it if both have a valid range)
		if (!Double.isNaN(m1) && !Double.isNaN(m2)) 
		{
			
			double p1 = (m1 - this.methMeans[dist]) / this.methSDs[dist];
			double p2 = (m2 - this.methMeans[dist]) / this.methSDs[dist];
			this.pearsonSums[dist] += (p1 * p2);
			this.counts[dist] += 2;
			
			// Mean and SD records. Count each seperately
			this.meanCounts[dist] ++;
			this.methTotals[dist] += m2;
			this.methSquareTotals[dist] += Math.pow(m2, 2);
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
	
	public double[] means()
	throws Exception
	{
		double[] countsCopy = Arrays.copyOf(this.meanCounts, this.meanCounts.length);
		// Replace zero counts with Nan, so they don't end up as Infinity
		for (int i = 0; i < countsCopy.length; i++)
		{
			if (this.meanCounts[i]<=0.0) countsCopy[i] = Double.NaN; 
		}
		double[] out = MatUtils.divVects(this.methTotals, countsCopy);
		return out;
		
//		double[] out = new double[this.counts.length];
//		double mean = this.methSummarizer.getValMean(false);
//		for (int i = 0; i < out.length; i++)
//		{
//			out[i] = mean;
//		}
//		return out;
	}

	public double[] stdevs()
	throws Exception
	{
		double[] means = this.means();
		double[] sds = new double[this.meanCounts.length];
		for (int i = 0; i < meanCounts.length; i++)
		{
			sds[i] = (this.meanCounts[i]==0) ? Double.NaN : Math.sqrt((this.methSquareTotals[i]/this.meanCounts[i])-Math.pow(means[i],2));
		}
		return sds;
		
//		double[] out = new double[this.counts.length];
//		double sd = this.methSummarizer.getValStdev();
//		for (int i = 0; i < out.length; i++)
//		{
//			out[i] = sd;
//		}
//		return out;
	}
	
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
			// Anything less than 20 can be unstable because we don't actually use the paired mean
			if (countsMinusOne[i]<=10.0) countsMinusOne[i] = Double.NaN; 
		}
		double[] pearsons = MatUtils.divVects(this.pearsonSums, countsMinusOne); 
		
		sb.append(ListUtils.excelLine(pearsons));
		//sb.append(ListUtils.excelLine(countsMinusOne));
		
		return sb.toString();
	}


	

}
