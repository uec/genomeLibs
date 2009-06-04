/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;
import java.util.zip.*;
import java.io.*;


import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;
import edu.usc.epigenome.genomeLibs.Counters.StringCounter;

/**
 * @author benb
 * 
 *
 */
public class APHandlerWigEmitter extends StringCounter implements AlignmentPosStreamHandler {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1641042070723916711L;
	
	
	// Types
	public static final int FIXED_STEP = 1;
	public static final int VARIABLE_STEP = 2;
	public static final int BEDGRAPH = 3;
	public static final int BED = 4;
	
	public final boolean strandedCounts = true;
	public final int defaultStrandedSpan = 0;
	
	protected int type = VARIABLE_STEP;
	protected String genome = null;
	protected int span = -1;
	protected String name = "trackName";
	protected String desc = "trackDesc";

	protected String lastChr = "none";
	protected boolean lastChrLegal = false;
	protected int lastCoord = -1;

	
	public PrintStream outstream = null;
	
	/**
	 *   Constructors
	 */
	public APHandlerWigEmitter(int inType, int inSpan, String inGenome, String inName, String inDesc) {
		try
		{
			initVars(inType, inSpan, inGenome, null, inName, inDesc);
		}
		catch (Exception e)
		{
			// Should actually never get here
			e.printStackTrace();
		}
	}

	public APHandlerWigEmitter(int inType, int inSpan, String inGenome, String outFilename,String inName, String inDesc)
	throws IOException 
	{
		initVars(inType, inSpan, inGenome, outFilename, inName, inDesc);
	}

	protected void initVars(int inType, int inSpan, String inGenome, String outFilename, String inName, String inDesc)
	throws IOException 
	{
		type = inType;
		genome = inGenome;
		name = inName;
		desc = inDesc;

		// If strandedCounts is set, the span doesn't apply because we are triangulating centerpoint
		span = (this.strandedCounts) ? this.defaultStrandedSpan : inSpan;
		
		if (outFilename == null)
		{
			outstream = System.out;
		}
		else
		{	
			OutputStream os = new FileOutputStream(outFilename);
			if (outFilename.toLowerCase().endsWith(".gz"))
			{
				os = new GZIPOutputStream(os);
			}

			outstream = new PrintStream(os);
		}
	}
	
	
	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	public void init() {
		outstream.printf("track type=wiggle_0 name=\"%s\" description=\"%s\"\n",name,desc);		
	}

	public void finish() {
    	System.err.println("Finishing emitter");
    	if (outstream != System.out) outstream.close();
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	public boolean streamElement(AlignmentPosStreamerPosition streamPos) 
	{
		boolean passes = true;
		double score = streamPos.getAvgScore(strandedCounts);
		String chr = streamPos.currentAp.getChr();
		int pos = streamPos.currentAp.getPos();
		
		try
		{
			switch (type)
			{
			case FIXED_STEP: streamFixedStep(chr, pos, score); break;
			case VARIABLE_STEP: streamVariableStep(chr, pos, score); break;
			//		case BED: streamFixedStep(chr, pos, score); break;
			//		case BEDGRAPH: streamFixedStep(chr, pos, score); break;
			default: System.err.println("Illegal APHandlerWigEmitter type: " + type); break;
			}
		}
		catch (Exception e)
		{
			System.err.println("APHandlerWigEmitter::streamElement() threw an exception: ");
			e.printStackTrace();
		}
		
		return passes;
	}
	
	protected void streamVariableStep(String chr, int pos, double score)
	throws Exception
	{
		
		// Filter out any non-USC chromosomes (like contam).
		boolean chrLegal = this.lastChrLegal;
		if (! (lastChr.equals(chr)))
		{
			// New chromosome
			// Is the chromosome legal
			chrLegal = GoldAssembly.chromExists(chr, this.genome);
			this.lastCoord = -1;
			
			if (chrLegal)
			{	
				outstream.printf("variableStep\tchrom=%s",chr);
				if (span > 1) outstream.printf("\tspan=%d", span);
				outstream.println();
			}
			
		}
			
		int coord = pos-Math.round(span/2);
		if (chrLegal)
		{
			if (coord <= 0)
			{
				coord = Math.max(1, lastCoord+1); // Sometimes we get multiple at 1
			}
			
			// UCSC also doesn't like any out of order coordinates (chrLegal also filters out the fake Ns at the end of the chrom)
			if (coord <= this.lastCoord)
			{
				throw new Exception("Why is " + chr + " coord " + this.lastCoord + " coming before coord " + pos);
			}

			// This one is unnecessary because we should always be *less* than the end of the chrom
			//coord = Math.min(coord, GoldAssembly.chromLengthStatic(chr, this.genome));
			
			outstream.printf("%d\t%f\n",coord, score);
		}	
		
		this.lastChr = chr;
		this.lastChrLegal = chrLegal;
		this.lastCoord = coord;
	}

	
	protected void streamFixedStep(String chr, int pos, double score)
	{
		
	}
	
	

}
