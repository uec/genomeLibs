package edu.usc.epigenome.testScripts;

import java.util.Iterator;
import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosIteratorMaqPileup;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosOptions;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APFilterCphs;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.DEPR__APHandlerWindowStatsCpGConcordance;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerWatsonThenCrick;

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
		AlignmentPosStreamer apStreamer = new AlignmentPosStreamerWatsonThenCrick(apIt, pre, post);
		
		// Add handlers, filters
		apStreamer.add(new APFilterCphs());
		DEPR__APHandlerWindowStatsCpGConcordance baseCounter = new DEPR__APHandlerWindowStatsCpGConcordance(100);
		apStreamer.add(baseCounter);
		
		// Run
		apStreamer.run();

	}
		
}
