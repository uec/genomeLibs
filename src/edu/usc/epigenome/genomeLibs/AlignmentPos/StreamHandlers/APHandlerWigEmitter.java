/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;
import java.util.zip.*;
import java.io.*;


import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;
import edu.usc.epigenome.genomeLibs.Counters.StringCounter;

/**
 * @author benb
 * 
 *
 */
public class APHandlerWigEmitter extends StringCounter implements AlignmentPosStreamHandler {

	// Types
	public static final int FIXED_STEP = 1;
	public static final int VARIABLE_STEP = 2;
	public static final int BEDGRAPH = 3;
	public static final int BED = 4;
	
	public int type = VARIABLE_STEP;
	public String genome = null;
	public int span = -1;
	public String name = "trackName";
	public String desc = "trackDesc";

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
		span = inSpan;
		genome = inGenome;
		name = inName;
		desc = inDesc;

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
		double score = streamPos.getAvgScore();
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
			// Is the chromosome legal
			chrLegal = GoldAssembly.chromExists(chr, this.genome);

			if (chrLegal)
			{	
				outstream.printf("variableStep\tchrom=%s",chr);
				if (span >= 0) outstream.printf("\tspan=%d", span);
				outstream.println();
			}
			
		}
		else
		{
			// UCSC also doesn't like any out of order coordinates
			if (pos <= this.lastCoord) throw new Exception("Why is " + chr + " coord " + this.lastCoord +
					" coming before coord " + pos);
		}

			
		if (chrLegal)
		{
			outstream.printf("%d\t%f\n",pos-Math.round(span/2), score);
		}	
		
		this.lastChr = chr;
		this.lastChrLegal = chrLegal;
		this.lastCoord = pos;
		
	}

	
	protected void streamFixedStep(String chr, int pos, double score)
	{
		
	}
	
	

}
