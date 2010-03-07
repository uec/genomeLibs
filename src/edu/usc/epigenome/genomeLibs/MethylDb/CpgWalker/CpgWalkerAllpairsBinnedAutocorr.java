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

public class CpgWalkerAllpairsBinnedAutocorr extends CpgWalkerAllpairs {

	protected double[] binEdges = null;
	protected TreeSet<Double> binTree = null;
	protected Map<Object, int[]> counters = null;
	
	
	public CpgWalkerAllpairsBinnedAutocorr(CpgWalkerParams inWalkParams,
			boolean inSamestrandOnly, double[] inBinEdges) {
		super(inWalkParams, inSamestrandOnly);
		
//		binEdges = inBinEdges;
		
		init(inBinEdges);
	}
	
	public CpgWalkerAllpairsBinnedAutocorr(CpgWalkerParams inWalkParams,
			boolean inSamestrandOnly, List<Double> inBinEdges) {
		super(inWalkParams, inSamestrandOnly);
		
//		binEdges = inBinEdges;
		
		double[] arr = new double[inBinEdges.size()];
		for (int i = 0; i < inBinEdges.size(); i++)
		{
			arr[i] = inBinEdges.get(i);
		}
		
		
		init(arr);
	}
	
	protected void init(double[] inBinEdges)
	{
		useSummarizers = false;
		
		// Initalize counters
		counters = new HashMap<Object,int[]>();
		int windSize = this.walkParams.maxScanningWindSize;

		// Setup bin tree
		binEdges = inBinEdges; // Save it for later
		binTree = new TreeSet<Double>();
		for (Double d : inBinEdges)
		{
			binTree.add(d);
		}
		
		for (int i = 0; i < (inBinEdges.length-1); i++)
		{
			double b1s = inBinEdges[i];
			//double b1e = (i==(inBinEdges.length-1)) ? b1s : inBinEdges[i+1];
			double b1e = inBinEdges[i+1];

			// singlet
			//System.err.printf("Initializing range <%f,%f>\n",b1s,b1e);
			counters.put(rangePairToKey(b1s,b1e,-1.0,-1.0), new int[windSize]);
			
			// And pairs (only up to the last one)
			for (int j = i; j < (inBinEdges.length-1); j++)
			{
				double b2s = inBinEdges[j];
//				double b2e = (j==(inBinEdges.length-1)) ? b2s : inBinEdges[j+1];
				double b2e = inBinEdges[j+1];
				
				System.err.printf("Initializing range <%f,%f,%f,%f>\n",b1s,b1e,b2s,b2e);
				counters.put(rangePairToKey(b1s,b1e,b2s,b2e), new int[windSize]);
			}
		}
	}

	@Override
	protected void recordPair(Cpg first, Cpg second)
	{
		// Get the distance
		int dist = second.chromPos - first.chromPos;
		

		Pair<Double,Double> r1 = getRange(first);
		Pair<Double,Double> r2 = getRange(second);
		
		int[] counter;
		
		// Add the singlet for each (only count it if both have a valid range)
		if ((r1 != null) && (r2 != null))
		{
			
			counter = this.counters.get(rangeSingletToKey(r1));
			counter[dist]++;

			counter = this.counters.get(rangeSingletToKey(r2));
			counter[dist]++;

			// And for the pair (we have to count it twice if we're in the same range)
			boolean sameRange = (r1.equals(r2));
			counter = this.counters.get(rangePairToKey(r1,r2));
			counter[dist]++;
			if (sameRange) counter[dist]++;
		}
		
	}
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#reset()
	 */
	@Override
	public void reset() {
		super.reset();
		init(this.binEdges);
	}

	
	@Override
	public String headerStr()
	{
		StringBuffer sb = new StringBuffer((int)1E5);
		
		sb.append("Range");
		for (int i = 0; i < this.walkParams.maxScanningWindSize; i++)
		{
			sb.append(",d");
			sb.append(i);
		}
		
		//System.err.println("Header="+sb.length());
		return sb.toString();
	}
	
	
	public String toCsvStrDumb()
	{
		StringBuffer sb = new StringBuffer((int)1E6);
		
		ListUtils.setDelim(",");
		for (Object key : counters.keySet())
		{
			sb.append(key.toString().replace(", ", "/"));
			sb.append(',');

			int[] counter = counters.get(key);
			sb.append(ListUtils.excelLine(counter));
			sb.append('\n');
		}
		
		
		return sb.toString();
	}
	
	@Override
	public String toCsvStr()
	throws Exception
	{
		StringBuffer sb = new StringBuffer((int)1E6);
		ListUtils.setDelim(",");

		for (int i = 0; i < (binEdges.length-1); i++)
		{
			double b1s = binEdges[i];
			double b1e = binEdges[i+1];
			double[] b1Counts = MatUtils.intVectToDouble(this.counters.get(rangePairToKey(b1s,b1e,-1.0,-1.0)));

			
			// And pairs (only up to the last one)
			for (int j = i; j < (binEdges.length-1); j++)
			{
				double b2s = binEdges[j];
				double b2e = binEdges[j+1];
				double[] b2Counts = MatUtils.intVectToDouble(this.counters.get(rangePairToKey(b2s,b2e,-1.0,-1.0)));
				
				double[] bothCounts = MatUtils.intVectToDouble(this.counters.get(rangePairToKey(b1s,b1e,b2s,b2e)));
						
				
				// Get it as a fraction of bin1
				sb.append(String.format("%.2f-%.2f->%.2f-%.2f",b1s,b1e,b2s,b2e,b1s,b1e));
				sb.append(',');
				sb.append(ListUtils.excelLine(MatUtils.divVects(bothCounts, b1Counts)));
				sb.append('\n');
				
				// And as a fraction of bin2 (unless they're identical)
				if (i!=j)
				{
				sb.append(String.format("%.2f-%.2f->%.2f-%.2f",b2s,b2e,b1s,b1e));
				sb.append(',');
				sb.append(ListUtils.excelLine(MatUtils.divVects(bothCounts, b2Counts)));
				sb.append('\n');
				}
				
			}
		}
		
		return sb.toString();
	}

	public Pair<Double,Double> getRange(Cpg cpg)
	{
		double val = cpg.fracMeth(true);
		if (Double.isNaN(val)) return null;
		
		return getRange(val);
	}
	
	public Pair<Double,Double> getRange(double d)
	{
		Double rs = binTree.floor(d);
		Double re = binTree.higher(d);
		
		// If it's the last one, we allow it to be equal
		if (re == null)
		{
			rs = binTree.lower(d);
			re = binTree.ceiling(d);
		}
		
		//System.err.printf("Val %.3f in bin <%.3f,%.3f>\n", d, rs, re);
		
		if ((rs == null) || (re == null))
		{
			System.err.printf("Error, couldn't get range for Val %0.3f\n\tbin <%0.3f,%0.3f>\n", d, rs, re);
			System.exit(1);
		}
		
		return Pair.from(rs, re);
	}
		
	static public Object rangePairToKey(Pair<Double,Double> r1, Pair<Double,Double> r2)
	{
		double r1s = Tuple.get1(r1);
		double r1e = Tuple.get2(r1);
		double r2s = Tuple.get1(r2);
		double r2e = Tuple.get2(r2);
		return rangePairToKey(r1s,r1e,r2s,r2e);
	}
	
	static public Object rangeSingletToKey(Pair<Double,Double> r1)
	{
		double r1s = Tuple.get1(r1);
		double r1e = Tuple.get2(r1);
		double r2s = -1.0;
		double r2e = -1.0;
		return rangePairToKey(r1s,r1e,r2s,r2e);
	}

	
	static public Object rangePairToKey(double r1s, double r1e, double r2s, double r2e)
	{
		// The smaller one always goes first 
		double lows, lowe, highs, highe;
		if (r1s <= r2s)
		{
			lows = r1s;
			lowe = r1e;
			highs = r2s;
			highe = r2e;
		}
		else
		{
			lows = r2s;
			lowe = r2e;
			highs = r1s;
			highe = r1e;
		}
		
		Object out = Quadruple.from(lows, lowe, highs, highe);
		return out;
	}
	

}
