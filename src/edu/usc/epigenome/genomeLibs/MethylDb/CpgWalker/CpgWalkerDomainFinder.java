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

	String curChr = null;
	PrintWriter pw = null;
	String lastChrom = "noChrom";
	
	List<GenomicRange> domains = new LinkedList<GenomicRange>();
	
	
	public CpgWalkerDomainFinder(CpgWalkerParams inWalkParams, String chr, PrintWriter inPw) {
		super(inWalkParams);
		pw = inPw;
	}

	
	/**
	 * @return the curChr
	 */
	public String getCurChr() {
		return curChr;
	}


	/**
	 * @param curChr the curChr to set
	 */
	public void setCurChr(String curChr) {
		this.curChr = curChr;
		super.newChrom();
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#reset()
	 */
	@Override
	public void reset() {
		super.reset();
	}


	@Override
	protected void processWindow()
	{
		boolean sameChrom = this.lastChrom.equalsIgnoreCase(this.getCurChr());

		
		// Does this window constitute a domain?
		int nCpgs = this.window.size();
		int s=-1, e=-1;

		if (nCpgs>0)
		{
			s = this.window.get(0).chromPos;
			e = this.window.get(nCpgs-1).chromPos;
		}
		
		
		boolean goodWind = ((nCpgs>0) && (nCpgs >= this.walkParams.minCpgs)); 
		goodWind &= this.windPasses();
		GenomicRange gr = new GenomicRange(this.getCurChr(), s, e);
		
		
		// Check if we overlap with the previous window
		boolean overlapsPrev = false;
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
			if (domains.size()>0)
			{
				GenomicRange lastGr = domains.get(domains.size()-1);
				String strandStr = (lastGr.getStrand() == StrandedFeature.NEGATIVE) ? "-" : "+";
				pw.append(MethylDbUtils.bedLine(lastGr.getChrom(), lastGr.getStart(), lastGr.getEnd(), strandStr, lastGr.getScore()));
				pw.println();
				domains.remove(0);
				dumpedLast = true;
			}
			
			// Then start new window if we're on a good one
			if (goodWind)
			{
				double score = this.windScore();
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
		
	
		this.lastChrom = this.curChr;
	}
	
	public List<GenomicRange> getDomains() {
		return domains;
	}


	public void setDomains(List<GenomicRange> domains) {
		this.domains = domains;
	}


	abstract protected boolean windPasses();
	abstract protected double windScore();
}
