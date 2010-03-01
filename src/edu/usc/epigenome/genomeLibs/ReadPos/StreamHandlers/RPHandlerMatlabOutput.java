package edu.usc.epigenome.genomeLibs.ReadPos.StreamHandlers;

import java.io.PrintStream;

import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;


/**
 * @author benb
 *
 * Outputs qual scores
 */
public class RPHandlerMatlabOutput implements ReadPosStreamHandler {

	PrintStream stream = null;
	int[][] qualMat = new int[5][1000];
	static final int EMPTY_VAL = Integer.MIN_VALUE;
	
	/**
	 * Streams to stdout
	 */
	public RPHandlerMatlabOutput() {
		stream = System.out;
	}

	/**
	 * Streams to a user specific stream
	 * 
	 * @param inStream
	 */
	public RPHandlerMatlabOutput(PrintStream inStream) {
		stream = inStream;
	}

	protected void emptyQualMat()
	{
		//System.err.println("Emptying int[" +qualMat.length+ "][" +qualMat[0].length+ "]");
		for (int i = 0; i < qualMat.length; i++)
		{
			for (int j = 0; j < qualMat[0].length; j++)
			{
				qualMat[i][j] = EMPTY_VAL;
			}
		}
	}
	
	protected void printQualMat(int maxCycle)
	{
		for (int i = 0; i < qualMat.length; i++)
		{
			for (int j = 0; j < maxCycle; j++)
			{
				if ((i>0)||(j>0)) stream.print(',');
				
				int val = qualMat[i][j];
				stream.print( (val==EMPTY_VAL) ? "nan" : val);
			}
		}
		stream.println(); // Newline
	}
	
	public void finish() {
		// TODO Auto-generated method stub

	}

	public void init() {
		// TODO Auto-generated method stub
		emptyQualMat();

	}

	public boolean streamElement(ReadPos currentRp) {

		int qual = currentRp.getQual();
		int cycle = currentRp.getCycle();
		int index = getSymbolIndex(currentRp);
		
		//System.err.println("\tIntermediate nuc of read, cycle " + currentRp.getCycle() + " (qual " + qual + ")");
	
		qualMat[index][cycle-1] = qual;
		
		if (currentRp.finalNucOfRead)
		{
			printQualMat(cycle);
			emptyQualMat();
		}

		
		return true;
	}
	
	/**
	 * @param rp
	 * @return a=0, c=1, t=2, g=3, n=4
	 */
	public static int getSymbolIndex(ReadPos rp)
	{
		int out = 4;
		switch (rp.getSymToken())
		{
		case 'A':
		case 'a':
			out = 0;
			break;
		case 'C':
		case 'c':
			out = 1;
			break;
		case 'G':
		case 'g':
			out = 2;
			break;
		case 'T':
		case 't':
			out = 3;
			break;
		case 'N':
		case 'n':
			out = 4;
			break;
		default:
			System.err.println("Don't understand symbol: " + rp.getSymToken());
			out = 4;
		break;		
		}
		return out;
	}

}
