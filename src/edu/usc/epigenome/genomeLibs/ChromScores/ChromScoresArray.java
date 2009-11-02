package edu.usc.epigenome.genomeLibs.ChromScores;

public class ChromScoresArray extends ChromScoresFast {

	public int downsamplingFactor = 1;
	
	
//	public ChromScoresArray() 
//	{
//		System.err.println("Initializing ChromScoresArray");
//		init();
//	}

	public ChromScoresArray(String genome) 
	{
		super(genome);
	}

	public ChromScoresArray(String genome, int inSmoothingFactor) 
	{
		super(genome);
		downsamplingFactor = inSmoothingFactor;
	}

	protected Object newChromArray(int chr_len)
	{
		//System.err.println("About to allocate array of " + chr_len + " floats");
		double[] out = new double[(int)Math.ceil((double)chr_len/(double)downsamplingFactor)+1];

		// Initialize to 0
		//for (int i = 0; i < chr_len; i++) out[i]=0.0;

		return out; 
	}
	

	protected Object addScoreToArray(Object array, int pos, Number score)
	{
		
		pos = Math.round((float)pos / (float)downsamplingFactor);
		
		if (pos >= ((double[])array).length)
		{
			System.err.println("Trying to add to a position ("+pos+") larger than the array size (" + ((double[])array).length + ")");
			System.exit(0);
		}
		else
		{
			double[] doubleArray = (double[])array;
			double incrementVal = (score.doubleValue() / (double)downsamplingFactor);
			doubleArray[pos] += incrementVal;
		}

		return array;
	}
	
	protected Number getArrayScore(Object array, int pos)
	{
		pos = Math.round(pos / downsamplingFactor);

		
		// Should we auto this to 0?
		double[] ar = (double[])array;
		Double out;
		if ((pos > 0) && (pos < ar.length))
		{
			out = new Double(ar[pos]);
		}
		else
		{
			out = new Double(0.0);
		}
			
		return out;
	}

	protected double[] getAllScores(Object array)
	{
		int s = minPos(array);
		int e = maxPos(array);
		int w = e - s + 1;
		
		double[] out = new double[w];
		System.arraycopy(array,s,out,0,w);
		//return (double[])array;
		return out;
	}

	 protected int minPos(Object array)
	 {
		 double[] ar = (double[])array;
		 
		 int ar_len = ar.length;
		 
		 return 0;
		 
//		 int first = 0;
//		 for ( first = 0; (first < ar_len) && (ar[first]<=0.0); first++)
//		 {
//		 }
//		 
//		 System.err.println("First = " + first);
//		 return first;
	 }

	
	 protected int maxPos(Object array)
	 {
		 double[] ar = (double[])array;
		 
		 int ar_len = ar.length;
		 System.err.println("ar_len = " + ar.length);
		 return ar.length - 1;
		 
//		 int last = 0;
//		 for ( last = (ar_len-1); (last > 0) && (ar[last]<=0.0); last--)
//		 {
//		 }
//		 
//		 System.err.println("Last = " + last);
//		 return last;
	 }

	

	
	
}
