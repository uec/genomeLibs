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
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APFilterFixedIntervals;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerGffEmitter;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerWigEmitter;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerWindowCounts;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer;

public class PileupToWig {

	// -c track cycles
	// -q track qual scores
	private static final String USAGE = "Usage: PileupToWig -output out.wig.gz -maxIdentical 1 -windSize 500 -stepSize -countEachBase file1.pileup file2.pileup ...";

	@Option(name="-genomeVers",usage="UCSC genome version code.  We need this to know the chromosome lengths.  Deafault \"hg18\"")
	private String genomeVers = "hg18";
	@Option(name="-output",usage="Output file. If suffix is \".gz\", it will be automatically gzip'd.  Deafault stdout")
	private String output = null;
	@Option(name="-name",usage="Track name.  Deafaults to filename")
	private String name = null;
	@Option(name="-desc",usage="Track description.  Deafaults to filename")
	private String desc = null;
    @Option(name="-windSize",usage="genomic window size (default 100)")
    private int windSize = 100;
    @Option(name="-type",usage="WIG format: FixedStep=1, VariableStep=2, BEDgraph=3, BED=4 (default 2)")
    private int type = APHandlerWigEmitter.VARIABLE_STEP;
    @Option(name="-stepSize",usage="genomic step size (default 100)")
    private int stepSize = 100;
    @Option(name="-maxIdentical",usage="Maximum reads with identical alignment positions (default infinite)")
    private int maxIdentical = 0;
    @Option(name="-countEachBase", usage="Count each base (default false, meaning count each read once)")
    private boolean countEachBase = false;
    
    // receives other command line parameters than options
    @Argument
    private List<String> arguments = new ArrayList<String>();

	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new PileupToWig().doMain(args);
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
        apos.trackPositions = true;
        apos.trackQuals = false;
        apos.maxIdentical = maxIdentical;
        apos.onlyFirstCycle = !countEachBase;

        APFilterFixedIntervals intervalFilter = new APFilterFixedIntervals(stepSize);
        {
        	
        	for (int i = 0; i < this.arguments.size(); i++)
        	{
        		String fn = (String)this.arguments.get(i);
            	String base = (new File(fn)).getName();
            	
            	String thisDesc = (desc != null) ? desc : base;

            	String thisName;
            	if (name == null)
            	{
            		thisName = base;
            	}
            	else if (this.arguments.size() > 1)
            	{
            		thisName = name + "." + base;
            	}
            	else 
            	{
            		thisName = name;
            	}
            	
            	APHandlerWigEmitter emitter = new APHandlerWigEmitter(type, windSize, genomeVers,output, thisName, thisDesc);

                	        		

        		// Create iterator, streamer
        		Iterator<AlignmentPos> apIt = new AlignmentPosIteratorMaqPileup(fn, apos);
        		int preWindSize = (int)Math.floor((double)windSize/2.0);
        		// Make the postWindow 1 shorter than pre because it does not include the current AP
        		int postWindSize = (int)Math.floor((double)(windSize-1)/2.0);

        		AlignmentPosStreamer apStreamer = new AlignmentPosStreamer(apIt,preWindSize, postWindSize);
        		apStreamer.add(intervalFilter);
        		apStreamer.add(emitter);

        		// Run
        		apStreamer.run();
        	}
        }
	}




		
}

