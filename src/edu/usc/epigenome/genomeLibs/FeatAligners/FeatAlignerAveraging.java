package edu.usc.epigenome.genomeLibs.FeatAligners;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;

import edu.usc.epigenome.genomeLibs.MatUtils;

public class FeatAlignerAveraging extends FeatAligner {

	// arr[0] fwTotalScores, arr[1] fwTotalFeats, arr[2] revTotalScores, arr[3] revTotalFeats
	double[][] arr;  
	
	/**
	 * @param flankSize
	 * @param zeroInit
	 */
	public FeatAlignerAveraging(int flankSize, boolean zeroInit) {
		super(flankSize, zeroInit);

		int nC = (flankSize*2) + 1;
		this.arr = new double[4][nC];
		
		MatUtils.initMat(arr, (zeroInit)  ? 0.0 : Double.NaN);
	}

	// Assumes that everything's already been flipped to the correct strand
	public void addAlignmentFast (double[] fwStrandScores, double[] revStrandScores)
	{
		for (int j = 0; j < fwStrandScores.length; j++)
		{
			double fwScore = fwStrandScores[j];
			double revScore = revStrandScores[j];
			this.addAlignmentPosFast(j, fwScore, revScore);
		}
	}
	
	// Assumes that everything's already been flipped to the correct strand
	protected void addAlignmentPosFast(int relPos, double fwScore, double revScore)
	{
		if (!Double.isNaN(fwScore))
		{
			this.arr[0][relPos]+= fwScore;
			this.arr[1][relPos]++;
		}

		if (!Double.isNaN(revScore))
		{
			this.arr[2][relPos]+= revScore;
			this.arr[3][relPos]++;
		}
	}

	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.FeatAligners.FeatAligner#toAverageFeatAligner()
	 */
	@Override
	public FeatAlignerAveraging toAverageFeatAligner() {
		return this;
	}


	@Override
	public void addAlignmentPos(int genomeRelPos, double fwStrandScore,
			double revStrandScore, String featName, String featChr,
			int featCoord, Strand featStrand) {
		int colInd = this.getColumnInd(genomeRelPos, featCoord, featStrand);
		this.addAlignmentPosFast(colInd,
				(featStrand == StrandedFeature.NEGATIVE) ? revStrandScore : fwStrandScore,
						(featStrand == StrandedFeature.NEGATIVE) ? fwStrandScore : revStrandScore);
	}


	
	
}
