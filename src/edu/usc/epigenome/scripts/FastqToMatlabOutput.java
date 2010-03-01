package edu.usc.epigenome.scripts;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.io.*;

import org.kohsuke.args4j.*;
import org.kohsuke.args4j.spi.*;

import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPosIteratorFastq;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPosOptions;
import edu.usc.epigenome.genomeLibs.ReadPos.StreamHandlers.RPHandlerMatlabOutput;
import edu.usc.epigenome.genomeLibs.ReadPos.StreamHandlers.RPHandlerReadPosCounts;
import edu.usc.epigenome.genomeLibs.ReadPos.StreamHandlers.ReadPosStreamHandler;
import edu.usc.epigenome.genomeLibs.ReadPos.Streamers.ReadPosStreamer;

public class FastqToMatlabOutput {

	// -c track cycles
	// -q track qual scores
	// -s use solexa qual
	private static final String USAGE = "Usage: FastqToMatlabOutput -illumina file1.fastq file2.fastq ...";
	
	
     @Option(name="-illumina",usage="Use solexa (64-scaled) quality scores, as opposed to Sanger default (33-scaled) (default false)")
    private boolean illumina = false;
    // receives other command line parameters than options
    @Argument
    private List<String> arguments = new ArrayList<String>();

	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new FastqToMatlabOutput().doMain(args);
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
		//rpos.minQualityScore = 0;
		rpos.trackPositions = true;
		rpos.trackQuals = true;
		rpos.positionQualsSolexaEncoding = this.illumina;
	    
		ReadPosStreamHandler matlabHandler = new RPHandlerMatlabOutput();
		for (int i = 0; i < this.arguments.size(); i++)
		{
			String fn = (String)this.arguments.get(i);

			// Create iterator, streamer
			Iterator<ReadPos> rpIt = new ReadPosIteratorFastq(fn, rpos);
			ReadPosStreamer rpStreamer = new ReadPosStreamer(rpIt);

			// Add handlers, filters
			rpStreamer.add(matlabHandler);

			// Run
			rpStreamer.run();
		}

	}
	



		
}

