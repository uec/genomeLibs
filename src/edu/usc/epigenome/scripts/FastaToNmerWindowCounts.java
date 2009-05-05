package edu.usc.epigenome.scripts;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.io.*;

import org.kohsuke.args4j.*;
import org.kohsuke.args4j.spi.*;

import edu.usc.epigenome.genomeLibs.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosIteratorFasta;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosIteratorMaqPileup;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosOptions;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APFilterNmer;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerWindowCounts;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer;

public class FastaToNmerWindowCounts {

	// -c track cycles
	// -q track qual scores
	private static final String USAGE = "Usage: FastaToNmerWindowCounts -nmer GAATG -strandSpecific -windSize 500 file1.fa file2.fa ...";
	
	
	@Option(name="-genomeVers",usage="UCSC genome version code.  We need this to know the chromosome lengths.  Deafault \"hg18\"")
	private String genomeVers = "hg18";
	@Option(name="-nmer",usage="nmer (default CG)")
    private String nmer = "CG";
    @Option(name="-windSize",usage="genomic window size (default 500)")
    private int windSize = 500;
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
    	new FastaToNmerWindowCounts().doMain(args);
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


		APHandlerWindowCounts counter = new APHandlerWindowCounts(windSize, strandSpecific, genomeVers);
		
		for (int i = 0; i < this.arguments.size(); i++)
		{
			String fn = (String)this.arguments.get(i);

			// Create iterator, streamer
			Iterator<AlignmentPos> apIt = new AlignmentPosIteratorFasta(fn);
			AlignmentPosStreamer apStreamer = new AlignmentPosStreamer(apIt, nmer.length(), nmer.length());
			
			apStreamer.add(new APFilterNmer(nmer, !strandSpecific));
			apStreamer.add(counter);

			// Run
			apStreamer.run();
		}	

		// Now output
		String description = "nmerCount," + nmer.toUpperCase() ;
		counter.excelOutput(description, System.out);
	}
	



		
}

