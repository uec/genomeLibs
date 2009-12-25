package edu.usc.epigenome.genomeLibs.FeatAligners;

import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;

import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;

public class FeatAlignerAveraging extends FeatAligner {

	// arr[0] fwTotalScores, arr[1] fwTotalFeats, arr[2] revTotalScores, arr[3] revTotalFeats
	double[][] arr;
	Set<GenomicRange> featsSeen;
	
	/**
	 * @param flankSize
	 * @param zeroInit
	 */
	public FeatAlignerAveraging(int flankSize, boolean zeroInit) {
		super(flankSize, zeroInit);

		int nC = (flankSize*2) + 1;
		this.arr = new double[4][nC];
		
		MatUtils.initMat(arr, (zeroInit)  ? 0.0 : Double.NaN);
		
		if (zeroInit)
		{
			// We want absolute numbers, so count the number of feats
			featsSeen = new HashSet<GenomicRange>();
		}
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

		if (this.zeroInit)
		{
			GenomicRange gr = new GenomicRange(featChr, featCoord, featCoord, featStrand);
			featsSeen.add(gr);
		}

	}


	public void matlabCsv(PrintWriter pw, boolean strandSpecific)
	{
		if (strandSpecific)
		{
			if (this.zeroInit)
			{
				// If it's zero init, it means we just want to divide by the
				// total.
				MatUtils.initMat(this.arr[1],(double)featsSeen.size());
				MatUtils.initMat(this.arr[3],(double)featsSeen.size());
			}
			MatUtils.matlabCsv(pw, this.arr);
		}
		else
		{
			
		}
	}

	
}
