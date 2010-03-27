package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import java.io.PrintWriter;
import java.util.*;

public class CpgWalkerPhasingFinder extends CpgWalker {
	
	protected boolean samestrandOnly = false;
	protected int period = 185;
	protected int halfPeriod = 92;
	protected PrintWriter pw = null;

	/**
	 * @param inWalkParams
	 */
	public CpgWalkerPhasingFinder(CpgWalkerParams inWalkParams, boolean inSamestrandOnly,
			int inPeriod, PrintWriter outWriter) {
		super(inWalkParams);
		this.samestrandOnly = inSamestrandOnly;
		this.period = inPeriod;
		this.halfPeriod = this.period / 2;
		System.err.printf("period=%d, half=%d\n",period, halfPeriod);
		this.pw = outWriter;
	}

	
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#newChrom()
	 */
	@Override
	protected void alertNewChrom() {
		super.alertNewChrom();
		//pw.printf("variableStep\tchrom=%s\n", curChr);
	}

	@Override
	protected void processWindow(List<Cpg> inWindow) {
		
		// Compare head node to each prior node in the window.
		int n = inWindow.size();
		if (n < 2) return;
		
		//System.err.println(this.windStr());
		
		Cpg head = inWindow.get(n-1);
		double total = 0.0;
		double numSeen = 0.0;
		for (int i = 0; i < (n-1); i++)
		{
			Cpg prior = inWindow.get(i);
			
			boolean process = true;
			if (samestrandOnly)
			{
				//System.err.println("Using same strand only");
				process = (prior.getStrand() == head.getStrand());
			}
			
			process &= ((head.chromPos-prior.chromPos)>10); // Don't count neighbors
			
			if (process)
			{
				//System.err.printf("\tProcessing pair: %d, %d\n", prior.chromPos, head.chromPos);
				
				if (prior.chromPos > head.chromPos)
				{
					// Only valid reason is if they're on different chromosomes.
					System.err.printf("Skipping out of order pair (new chrom?): %d, %d\n", prior.chromPos, head.chromPos);
					this.reset();
					return;
				}
				else if (prior.chromPos == head.chromPos)
				{
					System.err.printf("Why did we get the same CpG twice (%d)??\n",head.chromPos);  
					(new Exception()).printStackTrace();
					System.exit(1);
				}
				else
				{
					double pairVal = pairVal(prior, head);
					total+=pairVal;
					numSeen++;
				}
			}
		}
		
		// Output a wig line
		double avg = total / numSeen;
		if (numSeen>=1 && ((avg>0.5) || (avg<-0.5)))
		{
			//this.pw.printf("%d,total=%.2f\n", head.chromPos, avg);
			//this.pw.printf("%d\t%d\n",head.chromPos,Math.round(100.0*avg));
			this.pw.printf("%d,%d\n",head.chromPos,Math.round(100.0*avg));
		}
		
		// And the superclass
		super.processWindow(inWindow);
	}
	
	protected double pairVal(Cpg first, Cpg second)
	{
		double m1 = first.fracMeth(true);
		double m2 = second.fracMeth(true);
		int dist = second.chromPos-first.chromPos;
		double mult = multiplier(dist);
		double d = Math.abs(m1-m2);
		double score = (1.0-d) * mult;
		
		
		//pw.printf("\t%d - %d (dist %d) = %.2f (%.2f * %.2f)\n", first.chromPos, second.chromPos, dist,score, 1.0-d, mult);
		return score;
	}
	
	protected double multiplier(int dist)
	{
		int withinDist = dist % this.period;
		
		int distToMidpoint = Math.abs(withinDist - halfPeriod);
		double distToMidpointRel = Math.min(((double)distToMidpoint / (double)halfPeriod), 1.0);
		
		double mult = 0.0;
		if ((distToMidpointRel<0.15) || (distToMidpointRel>0.85))
		{
			mult = (2.0 * distToMidpointRel) - 1.0;
		}
		
		return mult;
	}
	
}
