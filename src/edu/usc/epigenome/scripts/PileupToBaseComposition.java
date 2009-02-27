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
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APFilterCpgs;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APFilterCphs;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerSymbolCounts;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerWatsonThenCrick;

public class PileupToBaseComposition {

	// -c track cycles
	// -q track qual scores
	private static final String USAGE = "Usage: PileupToBaseComposition -refComposition -additionalDesc DescritionTag -minQual 30 -cgonly -chonly -cycles -quals file1.pileup file2.pileup ...";
	
	
    @Option(name="-minQual",usage="minimum quality score (default 0)")
    private int minQual = 0;
    @Option(name="-cycles",usage="Store cycle info (default false)")
    private boolean cycles = false;
    @Option(name="-quals", usage="Store quality scores (default false)")
    private boolean quals = false;
    @Option(name="-cgonly", usage="Store quality scores (default false)")
    private boolean cgonly = false;
    @Option(name="-chonly", usage="Store quality scores (default false)")
    private boolean chonly = false;
    @Option(name="-refComposition",usage="Take the base composition of the REFERENCE sequence (useful for selection artifacts)")
    private boolean refComposition = false;
   @Option(name="-additionalDesc", usage="Additional description tag")
    private String additionalDesc = null;

    
    // receives other command line parameters than options
    @Argument
    private List<String> arguments = new ArrayList<String>();

	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new PileupToBaseComposition().doMain(args);
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
		apos.trackPositions = this.cycles;
		apos.trackQuals = this.quals;
		apos.maxIdentical = 0;

		APHandlerSymbolCounts counter = new APHandlerSymbolCounts(refComposition);
		
		for (int i = 0; i < this.arguments.size(); i++)
		{
			String fn = (String)this.arguments.get(i);

			// Create iterator, streamer
			Iterator<AlignmentPos> apIt = new AlignmentPosIteratorMaqPileup(fn, apos);
			AlignmentPosStreamer apStreamer = new AlignmentPosStreamerWatsonThenCrick(apIt, 1, 1);

			// Add handlers, filters
			if (chonly) apStreamer.add(new APFilterCphs());
			if (cgonly) apStreamer.add(new APFilterCpgs());
			apStreamer.add(counter);

			// Run
			apStreamer.run();

		}	

		// Now output
		String description = "postAlignment";
		if (additionalDesc != null) description += "_" + additionalDesc;
		if (chonly) description += ".CpH";
		if (cgonly) description += ".CpG";
		if (refComposition) description += ".RefComposition";
		
		System.out.print(counter.excelOutput(description));

	}
	



		
}

