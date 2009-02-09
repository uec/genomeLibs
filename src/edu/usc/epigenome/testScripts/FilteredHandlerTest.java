package edu.usc.epigenome.testScripts;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.io.*;

import org.kohsuke.args4j.*;
import org.kohsuke.args4j.spi.*;

import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.*;

public class FilteredHandlerTest {

	// -c track cycles
	// -q track qual scores
	private static final String USAGE = "Usage: PileupToDuplicateReport file1.pileup file2.pileup ...";
	
	
  
    // receives other command line parameters than options
    @Option(name="-maxIdentical",usage="Maximum reads with identical alignment positions (default infinite)")
    private int maxIdentical = 0;
    @Argument
    private List<String> arguments = new ArrayList<String>();

	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new FilteredHandlerTest().doMain(args);
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
		apos.maxIdentical = maxIdentical;
		apos.onlyFirstCycle = false;
		
		
		APHandlerDepthCounts cpgCounts =  new APHandlerDepthCounts();
		APHandlerDepthCounts cphCounts =  new APHandlerDepthCounts();
		APHandlerDepthCounts cCounts =  new APHandlerDepthCounts();

		for (int i = 0; i < this.arguments.size(); i++)
		{
			String fn = (String)this.arguments.get(i);

			// Create iterator, streamer
			Iterator<AlignmentPos> apIt = new AlignmentPosIteratorMaqPileup(fn, apos);
			AlignmentPosStreamer apStreamer = new AlignmentPosStreamer(apIt, 1, 1);
			apStreamer.add(new AlignmentPosStreamFilteredHandler(new APFilterCphs(),cphCounts));
			apStreamer.add(new AlignmentPosStreamFilteredHandler(new APFilterCpgs(),cpgCounts));
			apStreamer.add(new AlignmentPosStreamFilteredHandler(new APFilterCytosines(),cCounts));

			// Run
			apStreamer.run();
		}	

		// Now output
		System.out.print(cpgCounts.excelOutput("cpgReadCounts"));
		System.out.print(cphCounts.excelOutput("cphReadCounts"));
		System.out.print(cCounts.excelOutput("cytosineReadCounts"));

	}
	



		
}

