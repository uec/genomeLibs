package edu.usc.epigenome.genomeLibs;

import java.util.Iterator;

public class AlignmentPosStreamerWatsonThenCrick extends AlignmentPosStreamer {

	
	private static final long serialVersionUID = -975892714971917231L;
	
	protected int realPreWindSize;
	protected int realPostWindSize;
	
	

	/**
	 * @param inApIt
	 * @param inPreWindSize
	 * @param inPostWindSize
	 * 
	 * We use the larger of the two windows so that we can do forward and reverse.
	 */
	public AlignmentPosStreamerWatsonThenCrick(Iterator<AlignmentPos> inApIt,
			int inPreWindSize, int inPostWindSize) {
		super(inApIt, 
				Math.max(inPreWindSize,inPostWindSize), 
				Math.max(inPreWindSize,inPostWindSize));
		realPreWindSize = inPreWindSize;
		realPostWindSize = inPostWindSize;
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosStreamer#processAp(edu.usc.epigenome.genomeLibs.AlignmentPos[], edu.usc.epigenome.genomeLibs.AlignmentPos, edu.usc.epigenome.genomeLibs.AlignmentPos[])
	 */
	@Override
	protected boolean processAp(AlignmentPos[] fivePrimeAps,
			AlignmentPos currentAp, AlignmentPos[] threePrimeAps) {

		boolean out = true;
		for (int i = 0; i <= 1; i++)
		{
			boolean fw = (i==0);
			
			
			AlignmentPos[] preAps = new AlignmentPos[realPreWindSize];
			AlignmentPos[] postAps = new AlignmentPos[realPostWindSize];
			
			// First fill up the preAps.  Make sure to clone because we use
			// removeRevStrandReads , which is destructive
			for (int j = 0; j < realPreWindSize; j++)
			{
				if (fw) 
				{
					preAps[j] = fivePrimeAps[(fivePrimeAps.length - realPreWindSize) + j].clone(); 
				}
				else
				{
					preAps[j] = threePrimeAps[realPreWindSize - j - 1].flipped(); 
				}
				preAps[j].removeRevStrandReads();
			}
			
			// Now the postAps
			for (int j = 0; j < realPostWindSize; j++)
			{
				if (fw) 
				{
					postAps[j] = threePrimeAps[j].clone(); 
				}
				else
				{
					postAps[j] = fivePrimeAps[fivePrimeAps.length - j - 1].flipped(); 
				}
				postAps[j].removeRevStrandReads();
			}		

			// Flip the current one if necessary
			AlignmentPos currentApDirectional = (fw) ? currentAp.clone() : currentAp.flipped();
			currentApDirectional.removeRevStrandReads();
			
//			System.err.print("Streaming " + realPreWindSize + ", " + realPostWindSize + ":\t"); 
//			System.err.println(AlignmentPos.getRefTokens(preAps) + 
//					"," + currentApDirectional.getRefToken() + "," + AlignmentPos.getRefTokens(postAps));
			out &= super.processAp(preAps, currentApDirectional, postAps);
		}

		return out;
	}
	
	
	

}
