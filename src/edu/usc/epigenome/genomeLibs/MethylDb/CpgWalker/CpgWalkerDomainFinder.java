package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.GenomicRange.*;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;

abstract public class CpgWalkerDomainFinder extends CpgWalker {

	PrintWriter pw = null;
	List<GenomicRange> domains = new LinkedList<GenomicRange>();
	
	
	public CpgWalkerDomainFinder(CpgWalkerParams inWalkParams, String chr, PrintWriter inPw) {
		super(inWalkParams);
		pw = inPw;
	}

	

	/**
	 * @param curChr the curChr to set
	 */
	public void setCurChr(String curChr) {
		// This doesn't hurt if the last one was already processed
		boolean dumpedLast = this.outputAndRemoveCurWind();
		
		super.setCurChr(curChr);
	}

	/**
	 * Outputs the last domain
	 */
	public void finishChr() {
		boolean dumpedLast = this.outputAndRemoveCurWind();
		super.finishChr();
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#reset()
	 */
	@Override
	public void reset() {
		super.reset();
	}


	@Override
	protected void processWindow(List<Cpg> inWindow)
	{

		
		// Does this window constitute a domain?
		int nCpgs = inWindow.size();
		int s=-1, e=-1;

		if (nCpgs>0)
		{
			s = inWindow.get(0).chromPos;
			e = inWindow.get(nCpgs-1).chromPos;
		}

//		System.err.println("\tprocessing window: " + s +
//				"-" + e + "(" + 
//				 (1 + e - s) + ")");

		
		// This is where we call the actual filter
		boolean goodWind = ((nCpgs>0) && (nCpgs >= this.walkParams.minScanningWindCpgs)); 
		goodWind &= this.windPasses(inWindow);
		GenomicRange gr = new GenomicRange(this.getCurChr(), s, e);
		
		
		// Check if we overlap with the previous window
		boolean overlapsPrev = false;
		boolean sameChrom = !this.onNewChrom();
		if (sameChrom && (domains.size()>0))
		{
			// Do we actually overlap?
			GenomicRange last = domains.get(domains.size()-1);
			overlapsPrev = gr.overlaps(last);

			// Double check that we're actually past it (which we should be.
			if (!overlapsPrev)
			{
				if (last.getEnd()>=gr.getStart())
				{
					Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe(
							String.format("Why is previous non-overlapping window ahead of current one?? Quitting.\nPrev=%s, Cur=%s\n"));
					System.exit(1);
				}
			}
			
			// If we overlap, and it's a good window, we merge into the previous
			if (overlapsPrev  && goodWind)
			{
				last.setEnd(e);
			}
		}
		
		// If we don't overlap the previous domain, we can process it
		boolean dumpedLast = false;
		if (!overlapsPrev)
		{
			// If we have a prior one, dump it
			dumpedLast = this.outputAndRemoveCurWind();

			
			// Then start new window if we're on a good one
			if (goodWind)
			{
				double score = this.windScore(inWindow);
				gr.setScore(score);
				domains.add(gr);
			}
		}

		
		//// DEBUGGING
//		System.err.printf("Processed window %s-%s, goodWind=%s, dumpLast=%s, ndomains=%d\n", 
//				s,e,goodWind,dumpedLast, domains.size());
//		for (int i = 0; i < domains.size(); i++)
//		{
//			System.err.printf("\tdomain%d=%d-%d\n", i,domains.get(i).getStart(), domains.get(i).getEnd() );
//		}
////		if (goodWind)
////		{
////			for (Cpg c : this.window)
////			{
////				System.err.printf("\t%s\n", c.toStringExpanded());
////			}
////		}
		
	
		// And the superclass
		super.processWindow(inWindow);

	}
	
	/**
	 * @return true if we had a window to dump
	 */
	protected boolean outputAndRemoveCurWind()
	{
		
		// If we have a prior one, dump it
		boolean foundOne = false;
		if (domains.size()>0 )
		{
			GenomicRange lastGr = domains.get(domains.size()-1);
			int len = lastGr.getEnd()-lastGr.getStart()+1;
			if (len >= this.walkParams.minOutputWindSize)
			{
				if (this.walkParams.debug) System.err.printf("Output domain: %s\n", lastGr.toString());
				String strandStr = (lastGr.getStrand() == StrandedFeature.NEGATIVE) ? "-" : "+";
				pw.append(MethylDbUtils.bedLine(lastGr.getChrom(), lastGr.getStart(), lastGr.getEnd(), strandStr, lastGr.getScore()));
				pw.println();
			}
			domains.remove(0);
			foundOne = true;
		}
		return foundOne;
	}
	
	
	public List<GenomicRange> getDomains() {
		return domains;
	}


	public void setDomains(List<GenomicRange> domains) {
		this.domains = domains;
	}


	abstract protected boolean windPasses(List<Cpg> inWindow);
	abstract protected double windScore(List<Cpg>inWindow);
}
