/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;
import java.util.zip.*;
import java.io.*;


import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

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
	
	public String lastChr = "none";
	public int lastCoord = -1;
	public int type = VARIABLE_STEP;
	public int span = -1;
	
	public PrintStream outstream = null;
	
	/**
	 *   Constructors
	 */
	public APHandlerWigEmitter(int inType, int inSpan) {
		type = inType;
		span = inSpan;
		outstream = System.out;
	}

	public APHandlerWigEmitter(int inType, int inSpan, String outFilename)
	throws IOException 
	{
		type = inType;
		span = inSpan;

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
		outstream.printf("track type=wiggle_0 name=%s description=%s\n","wig1","descr");		
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
		
		switch (type)
		{
		case FIXED_STEP: streamFixedStep(chr, pos, score); break;
		case VARIABLE_STEP: streamVariableStep(chr, pos, score); break;
//		case BED: streamFixedStep(chr, pos, score); break;
//		case BEDGRAPH: streamFixedStep(chr, pos, score); break;
		default: System.err.println("Illegal APHandlerWigEmitter type: " + type); break;
		}
		
		return passes;
	}
	
	protected void streamVariableStep(String chr, int pos, double score)
	{
		if (! (lastChr.equals(chr)))
		{
			outstream.printf("variableStep\tchrom=%s",chr);
			if (span >= 0) outstream.printf("\tspan=%d", span);
			outstream.println();
		}
		
		outstream.printf("%d\t%f\n",pos-Math.round(span/2), score);
		
		this.lastChr = chr;
		this.lastCoord = pos;
	}

	
	protected void streamFixedStep(String chr, int pos, double score)
	{
		
	}
	
	

}
