/**
 * 
 */
package edu.usc.epigenome.genomeLibs.TrackFiles;

import java.io.*;
import java.nio.IntBuffer;

/**
 * @author benb
 *
 */
public class TrackFileRandomAccess extends TrackFile {

	protected RandomAccessFile raFile = null;
	
	/**
	 * @param myFile
	 * @param inGenome
	 * @param writable
	 */
	public TrackFileRandomAccess(File inFile, String inGenome, boolean writable, boolean initializeFullLength)
	throws IOException
	 {
		super(inFile, inGenome, writable, initializeFullLength);
		
		String mode = (writable) ? "rw" : "r";
		raFile = new RandomAccessFile(inFile, mode);
		
		//System.err.println("Genome size: " + this.getGenomeSize());

		if (writable)
		{
			raFile.setLength(this.getGenomeSize() * 4); // Because we use ints
		}
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TrackFiles.TrackFile#getValsRaw(java.lang.String, long, int, int[])
	 */
	@Override
	public void getValsRaw(String chrom, long startOffsetGlobal, int len,
			int[] out) throws IOException {

		long startFileOffset = startOffsetGlobal * 4; // Because it's an int
		raFile.seek(startFileOffset);

		System.err.println("Reading raw\t" + chrom + "\t" + startOffsetGlobal + "\t" + len);
		System.err.println("\tSeeking to " + startFileOffset);

		for (int i = 0; i < len; i++)
		{
			out[i] = raFile.readInt();
		}
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TrackFiles.TrackFile#writeChromRaw(java.lang.String, int[])
	 */
	@Override
	public void writeRaw(String chrom, long startOffsetGlobal, int[] vals)
	throws IOException
	{
		long startFileOffset = startOffsetGlobal * 4; // Because it's an int
		System.err.println("Writing raw\t" + chrom + "\t" + startOffsetGlobal + "\t" + vals.length);
		
		
		System.err.println("\tSeeking to " + startFileOffset);
		raFile.seek(startFileOffset);
		for (int i = 0; i < vals.length; i++)
		{
			if ((i%1E6)==0) System.err.println("\tOn line " + i);
			raFile.writeInt(vals[i]);
		}
		
		
	}

}
