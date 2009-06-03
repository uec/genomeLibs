package edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers;

import java.util.Iterator;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;

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
	protected boolean processAp(AlignmentPosStreamerPosition streamPos) {

		boolean out = true;
		for (int i = 0; i <= 1; i++)
		{
			boolean fw = (i==0);
			
			// Flip the current one if necessary
			AlignmentPos currentApDirectional = (fw) ? streamPos.currentAp.clone() : streamPos.currentAp.flipped();
			currentApDirectional.removeRevStrandReads();

//			// Only process if we have a non-zero depth after removing.
//			if (currentApDirectional.getTotalDepth() > 0)
//			{

				AlignmentPos[] preAps = new AlignmentPos[realPreWindSize];
				AlignmentPos[] postAps = new AlignmentPos[realPostWindSize];

				// First fill up the preAps.  Make sure to clone because we use
				// removeRevStrandReads , which is destructive
				for (int j = 0; j < realPreWindSize; j++)
				{
					if (fw) 
					{
						preAps[j] = streamPos.priorAps[(streamPos.priorAps.length - realPreWindSize) + j].clone(); 
						preAps[j].setStrand(StrandedFeature.POSITIVE);
					}
					else
					{
						preAps[j] = streamPos.nextAps[realPreWindSize - j - 1].flipped(); 
						preAps[j].setStrand(StrandedFeature.NEGATIVE);
					}
					preAps[j].removeRevStrandReads();
				}

				// Now the postAps
				for (int j = 0; j < realPostWindSize; j++)
				{
					if (fw) 
					{
						postAps[j] = streamPos.nextAps[j].clone(); 
						postAps[j].setStrand(StrandedFeature.POSITIVE);
					}
					else
					{
						postAps[j] = streamPos.priorAps[streamPos.priorAps.length - j - 1].flipped(); 
						postAps[j].setStrand(StrandedFeature.NEGATIVE);
					}
					postAps[j].removeRevStrandReads();
				}		


				//			System.err.print("Streaming " + realPreWindSize + ", " + realPostWindSize + ":\t"); 
				//			System.err.println(AlignmentPos.getRefTokens(preAps) + 
				//					"," + currentApDirectional.getRefToken() + "," + AlignmentPos.getRefTokens(postAps));
				AlignmentPosStreamerPosition newStreamPos = new AlignmentPosStreamerPosition(preWindSize, postWindSize);
				newStreamPos.priorAps = preAps;
				newStreamPos.currentAp = currentApDirectional;
				newStreamPos.nextAps = postAps;
				out &= super.processAp(newStreamPos);
			}
//		}

		return out;
	}
	
	
	

}
