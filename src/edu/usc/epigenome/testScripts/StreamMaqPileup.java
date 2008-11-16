package edu.usc.epigenome.testScripts;

import java.util.Iterator;
import edu.usc.epigenome.genomeLibs.*;

public class StreamMaqPileup {

	private static final String C_USAGE = "Usage: StreamMaqPileup max_identical pre post file1.pileup";
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception {

		if (args.length < 4)
		{
			System.err.println(C_USAGE);
			System.exit(0);
		}

		AlignmentPosOptions apos = new AlignmentPosOptions();
		String fn = null;
		int pre=0, post=0;
		try
		{
			apos.maxIdentical = Integer.parseInt(args[0]);
			pre = Integer.parseInt(args[1]);
			post = Integer.parseInt(args[2]);
			fn = args[3];
		}
		catch (Exception e)
		{
			System.err.println(C_USAGE + "\n" + e);
			System.exit(0);
		}
			
		// Create iterator, streamer
		Iterator<AlignmentPos> apIt = new AlignmentPosIteratorMaqPileup(fn, apos);
		AlignmentPosStreamer apStreamer = new AlignmentPosStreamer(apIt, pre, post);
		
		// Add handlers, filters
		APHandlerSymbolCounts baseCounter = new APHandlerSymbolCounts();
		apStreamer.add(baseCounter);
		
		// Run
		apStreamer.run();
		
		// Now output 
		System.out.print(baseCounter.excelOutput());

	}
		
}
