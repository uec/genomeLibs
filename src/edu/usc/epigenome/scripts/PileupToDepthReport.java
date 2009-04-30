package edu.usc.epigenome.scripts;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.io.*;

import org.kohsuke.args4j.*;
import org.kohsuke.args4j.spi.*;

import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerDepthCounts;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerDepthCountsRandomSubset;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer;

public class PileupToDepthReport {

	// -c track cycles
	// -q track qual scores
	private static final String USAGE = "Usage: PileupToDuplicateReport file1.pileup file2.pileup ...";
	
	
  
    // receives other command line parameters than options
    @Option(name="-maxIdentical",usage="Maximum reads with identical alignment positions (default infinite)")
    private int maxIdentical = 0;
    @Option(name="-randomSubset",usage="Use a ranomized subset of N read positions (N reads unless -countEachBase is set)")
    private int randomSubset = 0;
    @Option(name="-randomSubsetNumTrials",usage="Number of randomized trials to average (only valid when randomSubset is also set) (default 1)")
    private int randomSubsetNumTrials = 1;
    @Option(name="-countEachBase", usage="Count each base (default false, meaning count each read once)")
    private boolean countEachBase = false;
    @Argument
    private List<String> arguments = new ArrayList<String>();

	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new PileupToDepthReport().doMain(args);
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
        	
        	if ((randomSubset == 0) && (randomSubsetNumTrials != 1))
        	{
        		throw new CmdLineException("Can not set -randomSubsetNumTrials without setting -randomSubset");
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
		apos.maxIdentical = maxIdentical;
		apos.onlyFirstCycle = !countEachBase;
		
		
		APHandlerDepthCounts counter = (this.randomSubset > 0) ? (new APHandlerDepthCountsRandomSubset(randomSubsetNumTrials, randomSubset)) : 
			(new APHandlerDepthCounts());
		
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

