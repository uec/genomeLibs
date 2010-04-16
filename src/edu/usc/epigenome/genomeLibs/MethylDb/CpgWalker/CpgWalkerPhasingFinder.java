package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import java.io.PrintWriter;
import java.util.*;

import org.usckeck.genome.ChromFeatures;

public class CpgWalkerPhasingFinder extends CpgWalker {
	
	protected boolean matlabStyle = false;
	protected boolean samestrandOnly = false;
	protected int period = 185;
	protected int halfPeriod = 92;
	protected int flank = 0;
	protected PrintWriter pws[] = null;
	protected int chrInt = 0;
	
	protected double[] totals = new double[2];
	protected double[] numSeen = new double[2];

	/**
	 * @param inWalkParams
	 */
	public CpgWalkerPhasingFinder(CpgWalkerParams inWalkParams, boolean inSamestrandOnly,
			int inPeriod, PrintWriter outWriterA, PrintWriter outWriterB, boolean inMatlabStyle,
			int inFlank) {
		super(inWalkParams);
		this.samestrandOnly = inSamestrandOnly;
		this.matlabStyle = inMatlabStyle;
		this.period = inPeriod;
		this.halfPeriod = this.period / 2;
		this.flank = inFlank;
		System.err.printf("period=%d, half=%d\n",period, halfPeriod);
		this.pws = new PrintWriter[2];
		this.pws[0] = outWriterA;
		this.pws[1] = outWriterB;
	}

	
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker#newChrom()
	 */
	@Override
	protected void alertNewChrom() {
		super.alertNewChrom();
		String s = String.format("variableStep\tchrom=%s\tspan=1\n", curChr);
		if (!matlabStyle) pws[0].append(s);
		if (!matlabStyle) pws[1].append(s);
		
		chrInt = (new ChromFeatures()).chrom_from_public_str(curChr);
	}

	@Override
	protected void processWindow(List<Cpg> inWindow) {
		
		// Compare head node to each prior node in the window.
		int n = inWindow.size();
		if (n < 2) return;
		
		//System.err.println(this.windStr());
		
		Cpg head = inWindow.get(n-1);
		this.totals[0] = this.totals[1] = this.numSeen[0] = this.numSeen[1] = 0.0;
		for (int i = 0; i < (n-1); i++)
		{
			Cpg prior = inWindow.get(i);
			
			boolean process = true;
			if (samestrandOnly)
			{
				//System.err.println("Using same strand only");
				process = (prior.getStrand() == head.getStrand());
			}
			
//			process &= ((head.chromPos-prior.chromPos)>10); // Don't count neighbors
			
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
					//double pairVal = pairVal(prior, head);
					//total+=pairVal;
					//numSeen++;
					
					double m1 = head.fracMeth(true);
					double m2 = prior.fracMeth(true);
					double diff = Math.abs(m1-m2);
					int dist = head.chromPos-prior.chromPos;
					if ( (dist>=(halfPeriod-flank)) && (dist<=(halfPeriod+flank)))
					{
						this.numSeen[0]++;
						this.totals[0] += diff;
					}
					else if ( (dist>=(period-flank)) && (dist<=(period+flank)))
					{
						this.numSeen[1]++;
						this.totals[1] += diff;
					}
					
				}
			}
		}
		
		double val0 = (this.numSeen[0]==0.0) ? Double.NaN : (this.totals[0] / this.numSeen[0]);
		double val1 = (this.numSeen[1]==0.0) ? Double.NaN : (this.totals[1] / this.numSeen[1]);
		if (Double.isNaN(val0) && Double.isNaN(val1))
		{
			
		}
		else if (matlabStyle)
		{
			this.pws[0].printf("%d,%d,%.2f,%.2f\n", chrInt, head.chromPos,val0,head.getCpgWeight());
			this.pws[1].printf("%d,%d,%.2f,%.2f\n", chrInt, head.chromPos,val1,head.getCpgWeight());
		}
		else
		{
			// Wiggle format
			this.pws[0].printf("%d\t%d\n",head.chromPos,Math.round(100.0 * val0));
			this.pws[1].printf("%d\t%d\n",head.chromPos,Math.round(100.0 * val1));
		}
		
//		// Output a wig line
//		double avg = total / numSeen;
//		if (numSeen>=1 && ((avg>0.5) || (avg<-0.5)))
//		{
//			//this.pw.printf("%d,total=%.2f\n", head.chromPos, avg);
//			//this.pw.printf("%d\t%d\n",head.chromPos,Math.round(100.0*avg));
//			this.pw.printf("%d,%d\n",head.chromPos,Math.round(100.0*avg));
//		}
		
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
