package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

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
		
		boolean dumpLast = false;
		if (!goodWind)
		{
			if (domains.size()>0) dumpLast = true;
		}
		else
		{
			double thisScore = this.windScore();
			// Create a new domain
			GenomicRange gr = new GenomicRange(this.getCurChr(), s, e);
			gr.setScore(thisScore);
			
			// Do we need to merge it with our last domain, or start a new one?
			boolean merged = false;
			if (sameChrom && (domains.size()>0))
			{
				GenomicRange last = domains.get(domains.size()-1);
				if (gr.overlaps(last))
				{
					merged = true;
					last.setEnd(e);
				}
			}
			
			if (!merged)
			{
				domains.add(gr);

				// Print out the last domain and get rid of it
				if (domains.size()>1)
				{
					// Don't dump last if this is the first one
					dumpLast = true;
				}
			}
		}
		
//		System.err.printf("Processed window %s-%s, goodWind=%s, dumpLast=%s\n", 
//				s,e,goodWind,dumpLast);
//		if (goodWind)
//		{
//			for (Cpg c : this.window)
//			{
//				System.err.printf("\t%s\n", c.toStringExpanded());
//			}
//		}
		
		
		// Do we need to dump the last domain?
		if (dumpLast)
		{
			GenomicRange lastGr = domains.get(0);
			String strandStr = (lastGr.getStrand() == StrandedFeature.NEGATIVE) ? "-" : "+";
			pw.append(MethylDbUtils.bedLine(lastGr.getChrom(), lastGr.getStart(), lastGr.getEnd(), strandStr, lastGr.getScore()));
			pw.println();
			domains.remove(0);
		}

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
