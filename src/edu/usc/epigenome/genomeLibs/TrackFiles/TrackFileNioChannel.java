/**
 * 
 */
package edu.usc.epigenome.genomeLibs.TrackFiles;

import java.io.*;
import java.nio.*;
import java.nio.channels.FileChannel;

/**
 * @author benb
 *
 *	From Javadoc for channel.map()
 *     For most operating systems, mapping a file into memory is more expensive than
 *     reading or writing a few tens of kilobytes of data via the usual read and write
 *     methods. From the standpoint of performance it is generally only worth mapping
 *     relatively large files into memory.
 *
 */
public class TrackFileNioChannel extends TrackFile {

	protected ByteBuffer byteBuf = null;
	
	/**
	 * @param myFile
	 * @param inGenome
	 * @param writable
	 */
	public TrackFileNioChannel(File inFile, String inGenome, boolean writable, boolean initializeFullLength)
	throws IOException
	 {
		super(inFile, inGenome, writable, initializeFullLength);
		
		String mode = (writable) ? "rw" : "r";
		RandomAccessFile raFile = new RandomAccessFile(inFile, mode);
		
        FileChannel channel = raFile.getChannel();
        FileChannel.MapMode channelMode = (writable) ? FileChannel.MapMode.READ_WRITE : FileChannel.MapMode.READ_ONLY;

        byteBuf = channel.map(channelMode, 0, (int)channel.size());
		
		
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TrackFiles.TrackFile#getValsRaw(java.lang.String, long, int, int[])
	 */
	@Override
	public void getValsRaw(String chrom, long startOffsetGlobal, int len,
			int[] out) throws IOException {

		long startFileOffset = startOffsetGlobal * 4; // Because it's a long

//		raFile.seek(startFileOffset);
//		for (int i = 0; i < len; i++)
//		{
//			out[i] = raFile.readInt();
//		}
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.TrackFiles.TrackFile#writeChromRaw(java.lang.String, int[])
	 */
	@Override
	public void writeRaw(String chrom, long startOffsetGlobal, int[] vals)
	throws IOException
	{
		long startFileOffset = startOffsetGlobal * 4; // Because it's a long

//		raFile.seek(startFileOffset);
//		for (int i = 0; i < vals.length; i++)
//		{
//			raFile.writeInt(vals[i]);
//		}
		
		
	}

}
