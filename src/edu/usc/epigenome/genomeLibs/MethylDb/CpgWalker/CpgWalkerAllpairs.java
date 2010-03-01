package edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;

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
	protected void processWindow() {
		
		// Compare head node to each prior node in the window.
		int n = this.window.size();
		if (n < 2) return;
		
		//System.err.println(windStr());
		
		Cpg head = this.window.get(n-1);
		for (int i = 0; i < (n-1); i++)
		{
			Cpg prior = this.window.get(i);
			
			boolean process = true;
			if (samestrandOnly)
			{
				process = (prior.getStrand() == head.getStrand());
			}
			
			if (process)
			{
				//System.err.printf("\tProcessing pair: %d, %d\n", prior.chromPos, head.chromPos);
				recordPair(prior, head);
			}
		}
	}
	
	protected abstract void recordPair(Cpg first, Cpg second);
}
