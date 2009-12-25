package edu.usc.epigenome.genomeLibs.FeatAligners;

import java.io.PrintWriter;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;

public abstract class FeatAligner {

	public int flankSize;
	protected boolean zeroInit = false;


	
	/**
	 * @param flankSize
	 */
	/**
	 * @param flankSize
	 * @param zeroInit If true, initialize with 0s, otherwise initialize with NaNs.
	 */
	public FeatAligner(int flankSize, boolean zeroInit) {
		super();
		this.flankSize = flankSize;
		this.zeroInit = zeroInit;
	}

	
	
	/**
	 * @param fwStrandScores Use Double.NaN for unscored positions
	 * @param revStrandScores Use Double.NaN for unscored positions
	 * @param featName
	 * @param featChr
	 * @param featCoord
	 * @param strand If negative, scores will be about the center
	 */
	
	
	abstract public void addAlignmentPos (int genomeRelPos, double fwStrandScore, double revStrandScore,
			String featName, String featChr, int featCoord, StrandedFeature.Strand featStrand);
	
	abstract public FeatAlignerAveraging toAverageFeatAligner();
	abstract public void matlabCsv(PrintWriter pw, boolean strandSpecific);


	
	/*** Non-abstract ***/
	
	public void addAlignment (double[] fwStrandScores, double[] revStrandScores, 
			String featName, String featChr, int featCoord, StrandedFeature.Strand featStrand)
	{
		int nC = fwStrandScores.length;
		for (int i = 0; i < nC; i++)
		{
			int genomeRelPos = (featStrand == StrandedFeature.NEGATIVE) ? (featCoord + this.flankSize - i) : 
					(featCoord - this.flankSize + i);
			double fwScore = fwStrandScores[i];
			double revScore = (revStrandScores==null) ? Double.NaN : revStrandScores[i];
			this.addAlignmentPos(genomeRelPos, fwScore, revScore, featName, featChr, featCoord, featStrand);
		}
	}

	
	public void addAlignment (double[] fwStrandScores, 
			String featName, String featChr, int featCoord, StrandedFeature.Strand featStrand)
	{
		this.addAlignment(fwStrandScores, null, featName, featChr, featCoord, featStrand);
	}

	public void addAlignmentPos (int genomeRelPos, double fwStrandScore,
			String featName, String featChr, int featCoord, StrandedFeature.Strand featStrand)
	{
		this.addAlignmentPos(genomeRelPos, fwStrandScore, Double.NaN, featName, featChr, featCoord, featStrand);
	}
	
	/** Utility functions, takes care of strand reversal and nan vs. 0 stuff **/
	
	public int getColumnInd(int genomeRelPos, int featCoord, StrandedFeature.Strand featStrand,
			int rangeStart, int rangeEnd)
	{
		int out;
		if (featStrand == StrandedFeature.NEGATIVE)
		{
			out = rangeEnd - genomeRelPos;
		}
		else
		{
			out = genomeRelPos - rangeStart;
		}
		
		return out;				
	}

	public int getColumnInd(int genomeRelPos, int featCoord, StrandedFeature.Strand featStrand)
	{
		int rangeStart = featCoord - this.flankSize; 
		int rangeEnd = featCoord + this.flankSize; 
		
		return this.getColumnInd(genomeRelPos, featCoord, featStrand, rangeStart, rangeEnd);
	}
	

	
}
