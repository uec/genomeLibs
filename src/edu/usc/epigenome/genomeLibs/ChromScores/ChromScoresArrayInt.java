package edu.usc.epigenome.genomeLibs.ChromScores;


public class ChromScoresArrayInt extends ChromScoresFast {

	public ChromScoresArrayInt(String genome) {
		super(genome);
		init(genome);
	}

	
	protected Object newChromArray(int chr_len)
	{
		System.err.println("About to allocate array of " + chr_len + " floats");
		int[] out = new int[chr_len+1];

		// Initialize to 0
		//for (int i = 0; i < chr_len; i++) out[i]=0.0;

		return out; 
	}
	

	protected Object addScoreToArray(Object array, int pos, Number score)
	{
		((int[])array)[pos] += score.intValue();
		return array;
	}
	
	protected Number getArrayScore(Object array, int pos)
	{
		// Should we auto this to 0?
		return new Integer(((int[])array)[pos]);
	}

	protected double[] getAllScores(Object array)
	{
		int s = minPos(array);
		int e = maxPos(array);
		int w = e - s + 1;
		
		double[] out = new double[w];
		for (int i = s; i <= e; i++)
		{
			out[i-s] = (double)((int[])array)[i];
		}
		return out;
	}

	 protected int minPos(Object array)
	 {
		 int[] ar = (int[])array;
		 
		 int ar_len = ar.length;
		 int first = 0;
		 for ( first = 0; (first < ar_len) && (ar[first]==0); first++)
		 {
		 }
		 
		 System.err.println("First = " + first);
		 return first;
	 }

	
	 protected int maxPos(Object array)
	 {
		 int[] ar = (int[])array;
		 
		 int ar_len = ar.length;
		 int last = 0;
		 for ( last = (ar_len-1); (last > 0) && (ar[last]==0); last--)
		 {
		 }
		 
		 System.err.println("Last = " + last);
		 return last;
	 }

	

	
	
}
