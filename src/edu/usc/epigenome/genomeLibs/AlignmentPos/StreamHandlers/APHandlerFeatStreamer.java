/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.RangeLocation;

import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.AlignmentPos.*;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition;

/**
 * @author benb
 * 
 *
 */
 abstract public class APHandlerFeatStreamer implements AlignmentPosStreamHandler {
	 
	 public String gtfFilename = null;
	 public int windSize = 0;
	 
	 public long totalBaseCoverageCountFw = 0; // A counter of all bases 
	 public long totalBaseCoverageCountRev = 0; // A counter of all bases 
	 
	 public ChromFeatures cf = null;
	 public String cfCurChr = null;

	/**
	 * 
	 */
	public APHandlerFeatStreamer(String inGtfFilename, int inWindSize) {
		gtfFilename = inGtfFilename;
		windSize = inWindSize;
	}



	/*
	 * Overridden StreamHandler functions(non-Javadoc)
	 */
	
	public void init() {
//		super.init();

	}

	public void finish() {
//		super.finish();
	}


	// ********** OVERRIDE *************
	
	
	/**
	 * @param streamPos
	 * @param pair
	 * @param rec
	 * @param cpgRelativeOffset
	 * @return true if the Feat used the CpG
	 */
	abstract protected boolean streamFeat(AlignmentPosStreamerPosition streamPos, GFFRecord rec, int apRelativeOffset);
	

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.AlignmentPosStreamHandler#streamElement(edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamerPosition)
	 */
	public boolean streamElement(AlignmentPosStreamerPosition streamPos) 
	{
		AlignmentPos cur = streamPos.currentAp;
		String chr = cur.getChr();
		int chrInd = (new ChromFeatures()).chrom_from_public_str(chr);
		int curPos = cur.getPos();
		
		// Update total base counts
		this.totalBaseCoverageCountFw += cur.getDepth(true);
		this.totalBaseCoverageCountRev += cur.getDepth(false);
				
		// Are we on a new chromosome? 
		// checkChr returns false if ChromFeatures can't load the chromosome
		boolean chromOk = checkChr(chr); 
		boolean out = chromOk;
		if (chromOk)
		{
			// Build the target region.
			int regStart = cur.getPos() - windSize;
			int regEnd = cur.getPos() + windSize;
			Location regLoc = new RangeLocation(regStart, regEnd);
			GFFEntrySet es = cf.coord_filtered_features(chrInd, regLoc, false);
			

			Iterator it = es.lineIterator();
			out = false;
			while (it.hasNext())
			{
				GFFRecord rec = (GFFRecord)it.next();
//				System.err.println("\t" + regLoc + "\t" + GFFUtils.gffBetterString(rec));
				
				// Get the offset (BASED ON START POSITION)
				int offset = curPos - rec.getStart();
				if (rec.getStrand() == StrandedFeature.NEGATIVE) offset *= -1;

				boolean outThis = this.streamFeat(streamPos, rec, offset);
				out |= outThis; // Include it if any of the feats use it.
			}
		}
		
		return out;
	}
	
	
	// checkChr returns false if ChromFeatures can't load the chromosome
	protected boolean checkChr(String newChr)
	{
		
		boolean out =true; 
		
		if ((this.cfCurChr == null) || (!cfCurChr.equals(newChr)))
		{
			int newChrInd = (new ChromFeatures()).chrom_from_public_str(newChr);

			try
			{
				if (newChrInd <= 0) throw new Exception ("ChromFeatures doesn't understand chromosome " + newChr);
				
				// !!!!! This will complain if we use a non-hg18 chromosome !!!!
				cf = new ChromFeatures(gtfFilename, false);
				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				cf.c_first_chrom_num = newChrInd;	
				cf.c_last_chrom_num = newChrInd;
				System.err.println("Loading " + newChr + " from " + gtfFilename);
				cf.populate();
				System.err.println("\tLoaded chrom " + newChr + "\n" + cf.num_features() + " feats. Building tree");
				cf.preload_coord_filtered_features(newChrInd);
				System.err.println("\tBuilt tree for chrom " + newChr);
				cf.C_VERBOSITY = -3;
				
				// Set state
				cfCurChr = newChr;
			}
			catch (Exception e)
			{
				out = false;
				System.err.println(e.getMessage());
				e.printStackTrace();
			}
		}

		return out;
	}
	
	public long getTotalBaseCoverageCount()
	{
		return this.totalBaseCoverageCountFw + this.totalBaseCoverageCountRev;
	}

}
