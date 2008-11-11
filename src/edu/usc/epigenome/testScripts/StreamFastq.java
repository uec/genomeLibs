package edu.usc.epigenome.testScripts;

import java.util.Iterator;
import edu.usc.epigenome.genomeLibs.*;

public class StreamFastq {

	private static final String C_USAGE = "Usage: StreamFastq file1.fastq";
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception {

		if (args.length < 1)
		{
			System.err.println(C_USAGE);
			System.exit(0);
		}

		ReadPosOptions rpos = new ReadPosOptions();
		String fn = null;
		try
		{
			fn = args[0];
		}
		catch (Exception e)
		{
			System.err.println(C_USAGE + "\n" + e);
			System.exit(0);
		}
			
		// Create iterator, streamer
		Iterator<ReadPos> rpIt = new ReadPosIteratorFastq(fn, rpos);
		ReadPosStreamer rpStreamer = new ReadPosStreamer(rpIt);
		
		// Add handlers, filters
		RPHandlerSymbolCounts baseCounter = new RPHandlerSymbolCounts();
		rpStreamer.add(baseCounter);
		
		// Run
		rpStreamer.run();
		
		// Now output 
		System.out.print(baseCounter.excelOutput());

	}
		
}
