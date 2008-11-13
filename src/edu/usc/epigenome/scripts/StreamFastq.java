package edu.usc.epigenome.scripts;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.io.*;

import org.kohsuke.args4j.*;
import org.kohsuke.args4j.spi.*;

import edu.usc.epigenome.genomeLibs.*;

public class StreamFastq {

	// -c track cycles
	// -q track qual scores
	// -s use solexa qual
	private static final String USAGE = "Usage: StreamFastq -minQual 30 -c(ycles) -q(uals) -s(olexa) file1.fastq file2.fastq ...";
	
	
    @Option(name="-minQual",usage="minimum quality score (default 0)")
    private int minQual = 0;
    @Option(name="-cycles",usage="Store cycle info (default false)")
    private boolean cycles = false;
    @Option(name="-quals", usage="Store quality scores (default false)")
    private boolean quals = false;
    @Option(name="-solexa",usage="Use solexa (64-scaled) quality scores (default true)")
    private boolean solexa = true;
    // receives other command line parameters than options
    @Argument
    private List<String> arguments = new ArrayList<String>();

	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new StreamFastq().doMain(args);
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



	
		ReadPosOptions rpos = new ReadPosOptions();
		rpos.minQualityScore = this.minQual;
		rpos.trackPositions = this.cycles;
		rpos.trackQuals = this.quals;
		rpos.positionQualsSolexaEncoding = this.solexa;
	    
		String fn = (String)this.arguments.get(0);

			
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

