package edu.usc.epigenome.testScripts;

import java.util.Iterator;
import edu.usc.epigenome.genomeLibs.*;

public class StreamMaqPileup {

	private static final String C_USAGE = "Usage: ReadMaqPileup max_identical file1.pileup";
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception {

		if (args.length < 2)
		{
			System.err.println(C_USAGE);
			System.exit(0);
		}

		AlignmentPosOptions apos = new AlignmentPosOptions();
		String fn = null;
		try
		{
			apos.maxIdentical = Integer.parseInt(args[0]);
			fn = args[1];
		}
		catch (Exception e)
		{
			System.err.println(C_USAGE + "\n" + e);
			System.exit(0);
		}
			
		
		// Create iterator, streamer
		Iterator<AlignmentPos> apIt = new AlignmentPosIteratorMaqPileup(fn, apos);
		AlignmentPosStreamer apStreamer = new AlignmentPosStreamer(apIt, 2, 2);
		
		// Add handlers, filters
		APHandlerBaseCounts baseCounter = new APHandlerBaseCounts();
		baseCounter.MAX_CYCLES = 36;
		apStreamer.add(baseCounter);
		
		// Run
		apStreamer.run();
		
		// Now output 
		System.out.println(baseCounter.excelOutput());

	}
		
}
