package edu.usc.epigenome.genomeLibs;

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

public class PileupToCpgTemplate {

	// -c track cycles
	// -q track qual scores
	private static final String USAGE = "Usage: PileupToCpgXXXX file1.pileup file2.pileup ...";
	
	
  
    // receives other command line parameters than options
	@Option(name="-minQual",usage="minimum quality score (default 0)")
	protected int minQual = 0;
	@Option(name="-maxIdentical",usage="Maximum reads with identical alignment positions (default infinite)")
    protected int maxIdentical = 0;
    @Option(name="-minDepth",usage="minimum read depth (default 0)")
    protected int minDepth = 0;
    @Option(name="-minDepthEachStrand",usage="minimum depth applies to each strand")
    protected boolean minDepthEachStrand = false;
    @Option(name="-windSize",usage="window size, for CpG dinucs windows (default 1)")
    protected int windSize = 1;
    @Option(name="-CtTransitionFreq",usage="C->T transition rate (default 0.005)")
    protected double CtTransitionFreq = 0.005; // From Schmidt 2008, Li 2009 (unpublished)
    @Option(name="-cpgTrackFilename",usage="binary track file for CpG density")
    protected String cpgTrackFilename = null;
    @Option(name="-filterOutSnps",usage="Remove likely SNP CpGs")
    protected boolean filterOutSnps = false;
   @Argument
   protected List<String> arguments = new ArrayList<String>();

	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new PileupToCpgTemplate().doMain(args);
    }
    
    public void checkArgs()
    throws Exception
    {
    	if( arguments.isEmpty() )
    	{
    		throw new CmdLineException("Must supply at least one input file");
    	}
   	
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
        	checkArgs();
        	
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
		apos.minQualityScore = minQual;
		apos.trackPositions = true;
		apos.trackQuals = false;
		apos.trackBisulfiteConversion = true;
		apos.trackSnps = true;
		apos.maxIdentical = maxIdentical;
		apos.onlyFirstCycle = false;
		apos.CtTransitionFrequency = CtTransitionFreq;
		
		
		for (int i = 0; i < this.arguments.size(); i++)
		{
			String fn = (String)this.arguments.get(i);

			// Create iterator, streamer
			Iterator<AlignmentPos> apIt = new AlignmentPosIteratorMaqPileup(fn, apos);
			AlignmentPosStreamer apStreamer = new AlignmentPosStreamer(apIt, windSize, windSize);
			
			apStreamer.add(new APFilterCpgs());
			
			if (minDepth>0) apStreamer.add(new APFilterMinDepth(minDepth,minDepthEachStrand));
			if (filterOutSnps) apStreamer.add(new APHandlerCpgFilterNonSnpCpgs());
			
			this.addHandlers(apStreamer);

			// Run
			apStreamer.run();
		}	

	}
	

	protected void addHandlers(AlignmentPosStreamer apStreamer)
	{
		System.err.println("Can not call base class PileupToCpgTemplate::addHandlers() directly");
		System.exit(1);
	}

		
}

