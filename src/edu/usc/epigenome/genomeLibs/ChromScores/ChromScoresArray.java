package edu.usc.epigenome.genomeLibs.ChromScores;

public class ChromScoresArray extends ChromScoresFast {

//	public ChromScoresArray() 
//	{
//		System.err.println("Initializing ChromScoresArray");
//		init();
//	}

	public ChromScoresArray(String genome) 
	{
		super(genome);
	}

	protected Object newChromArray(int chr_len)
	{
		System.err.println("About to allocate array of " + chr_len + " floats");
		double[] out = new double[chr_len+1];

		// Initialize to 0
		//for (int i = 0; i < chr_len; i++) out[i]=0.0;

		return out; 
	}
	

	protected Object addScoreToArray(Object array, int pos, Number score)
	{
		
		if (pos >= ((double[])array).length)
		{
			System.err.println("Trying to add to a position ("+pos+") larger than the array size (" + ((double[])array).length + ")");
		}
		else
		{
			((double[])array)[pos] += score.doubleValue();
		}

		return array;
	}
	
	protected Number getArrayScore(Object array, int pos)
	{
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
		return out;
	}

	 protected int minPos(Object array)
	 {
		 double[] ar = (double[])array;
		 
		 int ar_len = ar.length;
		 int first = 0;
		 for ( first = 0; (first < ar_len) && (ar[first]<=0.0); first++)
		 {
		 }
		 
		 System.err.println("First = " + first);
		 return first;
	 }

	
	 protected int maxPos(Object array)
	 {
		 double[] ar = (double[])array;
		 
		 int ar_len = ar.length;
		 int last = 0;
		 for ( last = (ar_len-1); (last > 0) && (ar[last]<=0.0); last--)
		 {
		 }
		 
		 System.err.println("Last = " + last);
		 return last;
	 }

	

	
	
}