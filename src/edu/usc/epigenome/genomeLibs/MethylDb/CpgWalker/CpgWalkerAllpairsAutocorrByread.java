package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.HashMap;
//import java.util.List;
//import java.util.Map;
//import java.util.TreeSet;

//import com.mallardsoft.tuple.Pair;
//import com.mallardsoft.tuple.Quadruple;
//import com.mallardsoft.tuple.Tuple;

import java.util.Collection;

import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgRead;

public class CpgWalkerAllpairsAutocorrByread extends CpgWalkerAllpairs {

//	protected double[] binEdges = null;
//	protected TreeSet<Double> binTree = null;
//	protected Map<Object, int[]> counters = null;
	
	protected static final int SAME_READ = 1;
	protected static final int DIFFERENT_READ = 2;
	protected static final int ANY_READ = 3;
	
	protected int readType = ANY_READ;
	
	protected boolean useOnlyCG = true;
	
	protected int[] nMM;
	protected int[] nMU;
	protected int[] nUU; 
	protected int[] nM;
	protected int[] nU;
	
	
	public CpgWalkerAllpairsAutocorrByread(CpgWalkerParams inWalkParams, boolean inSamestrandOnly, boolean inOppstrandOnly, boolean inSameRead, boolean inDifferentRead) {
		super(inWalkParams, inSamestrandOnly, inOppstrandOnly);
		
		if (inSameRead && inDifferentRead)
		{
			System.err.println("CpgWalkerAllpairsAutocorrByread::constructor can not be called with both inSameRead and inDifferentRead set");
			(new Exception()).printStackTrace();
			System.exit(1);
		}
		else if (inSameRead)
		{
			readType = SAME_READ;
		}
		else if (inDifferentRead)
		{
			readType = DIFFERENT_READ;
		}
		else
		{
			readType = ANY_READ;
		}
		
		init();
	}
	


	/**
	 * This makes it reducible for map-reduce
	 * 
	 * @param cpgWalkerAllpairsAutocorrByread
	 * @param cpgWalkerAllpairsAutocorrByread2
	 * @return
	 */
	public static CpgWalkerAllpairsAutocorrByread merge(
			CpgWalkerAllpairsAutocorrByread a,
			CpgWalkerAllpairsAutocorrByread b) {
		
		System.err.printf("Merging Autcorr(%d) + Autocorr(%d)...", a.totalCount(), b.totalCount());
		
		CpgWalkerAllpairsAutocorrByread out = a;
		
		System.err.printf("Result has Autocorr(%d)\n", out.totalCount());
		return out;
	}

	
	
	protected void init()
	{
		// Initalize counters
		int windSize = this.walkParams.maxScanningWindSize;
		nMM = new int[windSize-1]; 
		nMU = new int[windSize-1]; 
		nUU = new int[windSize-1]; 
		nM = new int[windSize-1]; 
		nU = new int[windSize-1]; 
	}

	
	
	
	public boolean useOnlyCG() {
		return useOnlyCG;
	}



	public void useOnlyCG(boolean useOnlyCG) {
		this.useOnlyCG = useOnlyCG;
	}



	@Override
	protected void recordPair(Cpg a, Cpg b)
	{
		// Get the distance
		int dist = b.chromPos - a.chromPos - 1;
		
//		Map<Integer,CpgRead> aCgReads = a.getReads();
//		Map<Integer,CpgRead> bCgReads = b.getReads();
		
		Collection<CpgRead> aCgReads = a.getReads().values();
		Collection<CpgRead> bCgReads = b.getReads().values();

		//System.err.println("Comparing CpG (A): " + a.toString() + "\tdist=" + dist);

		for (CpgRead aCgRead : aCgReads)
		{
			boolean aCg = aCgRead.validCg( (this.walkParams.methylParams==null) || this.walkParams.methylParams.getUseNonconversionFilter());
			//if (!aCg) System.err.println("Got an uncounted cytosine (A): " + a.toString() + "\n\t" + aCgRead.toString()); // COMMENT OUT
			if (!useOnlyCG || aCg)
			{
				int aReadId = aCgRead.readId;
				boolean aMeth = aCgRead.meth(true);

				for (CpgRead bCgRead : bCgReads)
				{
					boolean bCg = bCgRead.validCg( (this.walkParams.methylParams==null) || this.walkParams.methylParams.getUseNonconversionFilter());
					//if (!bCg) System.err.println("Got uncounted cytosine (B): " + b.toString() + "\n\t" +  bCgRead.toString());  // COMMENT OUT
					if (!useOnlyCG || bCg)
					{
						int bReadId = bCgRead.readId;
						boolean identicalRead = (aReadId == bReadId);

						boolean include = (readType==ANY_READ) || (identicalRead&&(readType==SAME_READ)) || (!identicalRead&&(readType==DIFFERENT_READ));
						if (include)
						{
							boolean bMeth = bCgRead.meth(true);

							if (aMeth) { nM[dist]++; } else { nU[dist]++; }
							if (bMeth) { nM[dist]++; } else { nU[dist]++; }

							if (aMeth && bMeth)
							{
								nMM[dist]+=2;
							}
							else if (!aMeth && !bMeth)
							{
								nUU[dist]+=2;
							}
							else
							{
								nMU[dist]++;
							}
							

						}
					}
				}
			}
		}
		

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
	public String headerStr()
	{
		StringBuffer sb = new StringBuffer((int)1E5);
		
		//sb.append(String.format("readType = %d,\tsamestrand=%s\n", this.readType, (this.samestrandOnly)?"true":"false"  ));
		sb.append(String.format("%s,%s,%s,%s,%s,%s","Distance","nM","nU","nMM","nMU","nUU"));
		
		//System.err.println("Header="+sb.length());
		return sb.toString();
	}
	
	
	public String toCsvStrDumb()
	{
		StringBuffer sb = new StringBuffer((int)1E6);
		
		
//		ListUtils.setDelim(",");
//		for (Object key : counters.keySet())
//		{
//			sb.append(key.toString().replace(", ", "/"));
//			sb.append(',');
//
//			int[] counter = counters.get(key);
//			sb.append(ListUtils.excelLine(counter));
//			sb.append('\n');
//		}
		
		
		return sb.toString();
	}
	
	@Override
	public String toCsvStr()
	throws Exception
	{
		StringBuffer sb = new StringBuffer((int)1E6);
		//ListUtils.setDelim(",");

		for (int i = 0; i < nMM.length; i++)
		{
			sb.append(String.format("%d,%d,%d,%d,%d,%d\n",i+1,nM[i],nU[i],nMM[i],nMU[i],nUU[i]));
		}
		
		return sb.toString();
	}

	
	public int totalCount()
	{
		int[][] mat = new int[][]{nM, nU, nMM, nMU, nUU};
		return MatUtils.sumAll(mat);
	}


	

}
