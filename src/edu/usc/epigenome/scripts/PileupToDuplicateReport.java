package edu.usc.epigenome.scripts;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.io.*;

import org.kohsuke.args4j.*;
import org.kohsuke.args4j.spi.*;

import edu.usc.epigenome.genomeLibs.*;

public class PileupToDuplicateReport {

	// -c track cycles
	// -q track qual scores
	private static final String USAGE = "Usage: PileupToDuplicateReport file1.pileup file2.pileup ...";
	
	
  
    // receives other command line parameters than options
    @Argument
    private List<String> arguments = new ArrayList<String>();

	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new PileupToDuplicateReport().doMain(args);
    }
    
	public void doMain(String[] args)
	throws Exception
	{


		CmdLineParser parser = new CmdLineParser(this);
        // if you have a wider console, you could increase the value;
        // here 80 is also the default
        parser.setUsageWidth(80);
        try
		{
        	parser.parseArgument(args);

        	if( arguments.isEmpty() )
        	{
        		throw new CmdLineException("Must supply at least one input file");
        	}
		}
        catch (CmdLineException e)
        {
        	System.err.println(e.getMessage());
            System.err.println(USAGE);
            // print the list of available options
            parser.printUsage(System.err);
            System.err.println();
            return;
        }



	
		AlignmentPosOptions apos = new AlignmentPosOptions();
		apos.minQualityScore = 0;
		apos.trackPositions = true;
		apos.trackQuals = false;
		apos.maxIdentical = 0;
		apos.onlyFirstCycle = true; // Track at the read level

		APHandlerDepthCounts counter = new APHandlerDepthCounts();
		
		for (int i = 0; i < this.arguments.size(); i++)
		{
			String fn = (String)this.arguments.get(i);

			// Create iterator, streamer
			Iterator<AlignmentPos> apIt = new AlignmentPosIteratorMaqPileup(fn, apos);
			AlignmentPosStreamer apStreamer = new AlignmentPosStreamer(apIt, 0, 0);
			apStreamer.add(counter);

			// Run
			apStreamer.run();
		}	

		// Now output
		String description = "identicalReadCounts";
		
		System.out.print(counter.excelOutput(description));

	}
	



		
}

