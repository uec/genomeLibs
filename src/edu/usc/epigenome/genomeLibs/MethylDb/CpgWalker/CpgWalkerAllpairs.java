package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import java.util.*;

public abstract class CpgWalkerAllpairs extends CpgWalker {
	
	protected boolean samestrandOnly = false;

	/**
	 * @param inWalkParams
	 */
	public CpgWalkerAllpairs(CpgWalkerParams inWalkParams, boolean inSamestrandOnly) {
		super(inWalkParams);
		this.samestrandOnly = inSamestrandOnly;
	}

	@Override
	protected void processWindow(List<Cpg> inWindow) {
		
		// Compare head node to each prior node in the window.
		int n = inWindow.size();
		if (n < 2) return;
		
		//System.err.println(windStr());
		
		Cpg head = inWindow.get(n-1);
		for (int i = 0; i < (n-1); i++)
		{
			Cpg prior = inWindow.get(i);
			
			boolean process = true;
			if (samestrandOnly)
			{
				System.err.println("Using same strand only");
				process = (prior.getStrand() == head.getStrand());
			}
			
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
				else
				{
					recordPair(prior, head);
				}
			}
		}
	}
	
	protected abstract void recordPair(Cpg first, Cpg second);
}
