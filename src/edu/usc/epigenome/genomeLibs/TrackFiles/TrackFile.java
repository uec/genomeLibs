/**
 * 
 */
package edu.usc.epigenome.genomeLibs.TrackFiles;

import java.io.*;

import org.apache.commons.collections.ListUtils;
import org.apache.commons.math.stat.StatUtils;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.MatUtils;

/**
 * @author benb
 *
 */
public abstract class TrackFile {

	protected File myFile = null;
	protected String genome = null;
	
	
	
	/**
	 * @param myFile
	 */
	public TrackFile(File myFile, String inGenome, boolean writable, boolean initializeFullLength)
	throws IOException
	{
		super();
		this.myFile = myFile;
		this.setGenome(inGenome);
	}
	
	
	public void setGenome(String inGenome)
	{
		this.genome = inGenome;
	}

	public String getGenome()
	{
		return this.genome;
	}
	
	public long getGenomeSize()
	{
		long out = GoldAssembly.getGenomeSize(this.getGenome());
		return out;
	}
	

	public int getVal(String chrom, int pos)
	throws IOException
	{
		int[] ar = getVals(chrom, pos, 1);
		return ar[0];
	}
	
	
	public double getValTotal(String chrom, int startPos, int len)
	throws IOException
	{
		System.err.println("Getting " + len + " vals");
		int[] all = getVals(chrom,startPos,len);

		System.err.println("Converting to double[]");
		double[] allD = edu.usc.epigenome.genomeLibs.ListUtils.intArrToDouble(all);

		System.err.println("Taking mean");
		double sum = StatUtils.sum(allD);
		
//		System.err.println(edu.usc.epigenome.genomeLibs.ListUtils.excelLine(all));

		return sum;
	}
	
	public double getValAvg(String chrom, int startPos, int len)
	throws IOException
	{
		System.err.println("Getting " + len + " vals");
		int[] all = getVals(chrom,startPos,len);

		System.err.println("Converting to double[]");
		double[] allD = edu.usc.epigenome.genomeLibs.ListUtils.intArrToDouble(all);

		System.err.println("Taking mean");
		double mean = StatUtils.mean(allD);
		
//		System.err.println(edu.usc.epigenome.genomeLibs.ListUtils.excelLine(all));

		return mean;
	}

	
	public int[] getVals(String chrom, int startPos, int len)
	throws IOException
	{
		int[] out = new int[len];
		this.getVals(chrom,startPos,len,out);
		return out;
	}
	
	public void getVals(String chrom, int startPos, int len, int[] out)
	throws IOException
	{
		long startOffset = GoldAssembly.getGlobalOffset(chrom, this.getGenome(), startPos);
		
		getValsRaw(chrom, startOffset, len, out);
	}

	public void writeChrom(String chrom, int[] vals)
	throws Exception
	{
		// Check that the input array has the right length.
		int expectedLen = GoldAssembly.chromLengthStatic(chrom, this.getGenome());
		
		if (expectedLen != vals.length)
		{
			throw new Exception("TrackFile::writeChrom() error: Chrom " + chrom + 
					" has expected length of " + expectedLen + " but got length " + vals.length);
		}
		
		long startOffset = GoldAssembly.getGlobalOffset(chrom, this.getGenome(), 0);
		this.writeRaw(chrom, startOffset, vals);
	}

	
	/*
	 * Abstract methods
	 */
	
	abstract public void getValsRaw(String chrom, long startOffsetGlobal, int len, int[] out) throws IOException;
	abstract public void writeRaw(String chrom, long startOffsetGlobal, int[] vals) throws IOException;
	
}
