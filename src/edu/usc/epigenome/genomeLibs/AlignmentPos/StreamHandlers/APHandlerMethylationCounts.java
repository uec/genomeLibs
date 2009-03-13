/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.*;

import edu.usc.epigenome.genomeLibs.GFFUtils;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosSnps;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosSnpsBisulfiteConverted;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;
import edu.usc.epigenome.genomeLibs.Counters.StringCounter;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;

/**
 * @author benb
 * 
 *
 */
public class APHandlerMethylationCounts extends StringCounter implements AlignmentPosStreamHandler {

	 public int fChgsStreamed = 0;
	 public int fChgsPassed = 0;
	 public int fChgReadsStreamed = 0;
	 public int fChgReadsPassed = 0;
	 protected boolean fReportCounts = true;
	 
	/**
	 * Constructor
	 */
	public APHandlerMethylationCounts() {
	}

	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	public void init() {
	}

	public void finish() {
		if (fReportCounts)
		{
		System.err.println("Num unique CHGs streamed: " + fChgsStreamed);
		System.err.println("Num unique CHGs streamed, passed: " + fChgsPassed);
		System.err.println("Num CHG reads streamed: " + fChgReadsStreamed);
		System.err.println("Num CHG reads streamed, passed: " + fChgReadsPassed);
		}
		
		System.out.println(this.excelOutput());
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamHandler#streamElement(java.util.LinkedList, edu.usc.epigenome.genomeLibs.AlignmentPos, java.util.LinkedList)
	 */
	public boolean streamElement(AlignmentPosStreamerPosition streamPos) 
	{
		boolean passes = true;
		
		AlignmentPosSnpsBisulfiteConverted ap = (AlignmentPosSnpsBisulfiteConverted)streamPos.currentAp;

		double meth = ap.getMethylatedFrac();
		if (Double.isNaN(meth))
		{
			this.increment("NaN");
		}
		else
		{
			double rounded = Math.floor(meth * 10) / 10;
			this.increment(Double.toString(rounded));
		}
		
		//System.err.println(GFFUtils.gffCsvLine(ap.toGff(true)));
		
		int depth = ap.getDepth(true);
		 fChgsStreamed++;
		 fChgReadsStreamed += depth;
		 if (passes)
		 {
			 fChgsPassed++;
			 fChgReadsPassed += depth;
		 }
		
		return passes;
	}

	
	
	

}
