package edu.usc.epigenome.scripts;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.io.*;

import org.kohsuke.args4j.*;
import org.kohsuke.args4j.spi.*;

import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosIteratorMaqPileup;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosOptions;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerWindowCounts;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer;

public class PileupToReadDepthWindows {

	// -c track cycles
	// -q track qual scores
	private static final String USAGE = "Usage: PileupToReadDepthWindows -strandSpecific -maxIdentical 1 -minQual 30 -windSize 500 -countEachBase file1.pileup file2.pileup ...";
	
	
    @Option(name="-minQual",usage="minimum quality score (default 0)")
    private int minQual = 0;
    @Option(name="-windSize",usage="genomic window size (default 500)")
    private int windSize = 500;
    @Option(name="-maxIdentical",usage="Maximum reads with identical alignment positions (default infinite)")
    private int maxIdentical = 0;
    @Option(name="-countEachBase", usage="Count each base (default false, meaning count each read once)")
    private boolean countEachBase = false;
    @Option(name="-strandSpecific", usage="Count reads on opposite strands separately (field 4 of output lists strand)")
    private boolean strandSpecific = false;
    
    // receives other command line parameters than options
    @Argument
    private List<String> arguments = new ArrayList<String>();

	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new PileupToReadDepthWindows().doMain(args);
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
		apos.minQualityScore = this.minQual;
		apos.trackPositions = true;
		apos.trackQuals = false;
		apos.maxIdentical = maxIdentical;
		apos.onlyFirstCycle = !countEachBase;

		APHandlerWindowCounts counter = new APHandlerWindowCounts(windSize, strandSpecific);
		
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
		String description = "readDepthWind";
		counter.excelOutput(description, System.out);
	}
	



		
}

