package edu.usc.epigenome.genomeLibs.Counters;

public class SNPByCycleCounter {

	static final int BUFFERLEN = 1000000;
	static final int DEFAULT_CYCLES = 130;

	int[][][] mMat = null;
	int mMaxCycle = 0;
	
	public SNPByCycleCounter() {
		mMat = new int[6][6][DEFAULT_CYCLES];
	}

	public SNPByCycleCounter(int numCycles) {
		mMat = new int[6][6][numCycles];
	}

	public void increment(char nucA, char nucB, int cycle)
	throws Exception
	{
		int indexA = nucToIndex(nucA);
		int indexB = nucToIndex(nucB);
		mMat[indexA][indexB][cycle]++;
		if (cycle>mMaxCycle) mMaxCycle = cycle;
	}
	
	public void increment(char nucA, char nucB, int cycle, int count)
	throws Exception
	{
		int indexA = nucToIndex(nucA);
		int indexB = nucToIndex(nucB);
		mMat[indexA][indexB][cycle] += count;
		if (cycle>mMaxCycle) mMaxCycle = cycle;
	}

	

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder(BUFFERLEN);
		
		int nN = mMat.length; // Number of nucleotides
		int nR = mMaxCycle;
		
		for (int r = 0; r < nR; r++)
		{
			for (int a = 0; a < nN; a++)
			{
				for (int b = 0; b < nN; b++)
				{
					if (a>0 || b>0) sb.append(',');
					sb.append( mMat[a][b][r] );
				}
			}
			sb.append('\n');
		}
		
		return sb.toString();
	}

	public static int nucToIndex(char nuc)
	throws Exception
	{
		int out = -1;
		
		switch (nuc)
		{
		case 'A': case 'a': out = 0; break;
		case 'C': case 'c': out = 1; break;
		case 'G': case 'g': out = 2; break;
		case 'T': case 't': out = 3; break;
		case 'N': case 'n': out = 4; break;
		case '-': out = 5; break;
		default: throw new Exception("Found an illegal nucleotide character: " + nuc);
		}

		return out;
	}
}
