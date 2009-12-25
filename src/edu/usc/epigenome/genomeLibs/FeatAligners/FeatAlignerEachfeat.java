package edu.usc.epigenome.genomeLibs.FeatAligners;

import java.io.PrintWriter;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;

import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;

public class FeatAlignerEachfeat extends FeatAligner {

	// i = type: arr[0] fwTotalScores, arr[1] revTotalScores
	// j = featNum: arr[0][5] = fwTotalScores for feat 6
	// k = coordinate: arr[0][5][350] = coord (350 - flank) relative to feature center.
	
	double[][][] arr;
	String[] featNames;
	GenomicRange[] featCoords;
	Map<GenomicRange,Integer> featinds = new TreeMap<GenomicRange,Integer>();
	int nFeatsSeen = 0;
	
	/**
	 * @param flankSize
	 * @param zeroInit
	 */
	public FeatAlignerEachfeat(int flankSize, boolean zeroInit, int nFeats) {
		super(flankSize, zeroInit);

		int nC = (flankSize*2) + 1;
		this.arr = new double[4][nFeats][nC];
		this.featCoords = new GenomicRange[nFeats];
		this.featNames = new String[nFeats];
		
		MatUtils.initMat(arr, (zeroInit)  ? 0.0 : Double.NaN);
	}

	
	public int numFeats()
	{
		return nFeatsSeen;
	}
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.FeatAligners.FeatAligner#addAlignmentPos(double, double, java.lang.String, java.lang.String, int, org.biojava.bio.seq.StrandedFeature.Strand)
	 */
	@Override
	public void addAlignmentPos(int genomeRelPos, double fwStrandScore, double revStrandScore,
			String featName, String featChr, int featCoord, Strand featStrand) {

		GenomicRange gr = new GenomicRange(featChr, featCoord, featCoord, featStrand);
		int featInd = this.getInd(gr, featName);
		int colInd = this.getColumnInd(genomeRelPos, featCoord, featStrand);

		// Flip the strand of the scores if features are flipped
		arr[0][featInd][colInd] = (featStrand==StrandedFeature.NEGATIVE) ? revStrandScore : fwStrandScore;
		arr[1][featInd][colInd] = (featStrand==StrandedFeature.NEGATIVE) ? fwStrandScore : revStrandScore;
	}




	@Override
	public String htmlChart(boolean strandSpecific) throws Exception{

		FeatAlignerAveraging av = this.toAverageFeatAligner();
		return av.htmlChart(strandSpecific);
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.FeatAligners.FeatAligner#toAverageFeatAligner()
	 */
	@Override
	public FeatAlignerAveraging toAverageFeatAligner() {

		FeatAlignerAveraging out = new FeatAlignerAveraging(this.flankSize, this.zeroInit);
		
		for (int j = 0; j < this.numFeats(); j++)
		{
			double[] fwScores = arr[0][j];
			double[] revScores = arr[1][j];
			out.addAlignmentFast(fwScores, revScores);
		}
		
		return out;
	}

	
	private int getInd(GenomicRange gr, String featName) {
		
		Integer ind = featinds.get(gr);
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(
				String.format("Adding new feature at pos=%d to ind=%d\n",gr.getStart(),ind));
		if (ind == null)
		{
			ind = new Integer(nFeatsSeen++);
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(
					String.format("Adding new feature at pos=%d to ind=%d\n",gr.getStart(),ind));
			featinds.put(gr, ind);
			this.featNames[ind.intValue()] = featName;
			this.featCoords[ind] = gr;
		}
		return ind;
	}

	public void matlabCsv(PrintWriter pw, boolean strandSpecific)
	{
		if (strandSpecific)
		{
			MatUtils.matlabCsv(pw, this.arr[0], this.numFeats(), 0);
			MatUtils.matlabCsv(pw, this.arr[1], this.numFeats(), 0);
		}
		else
		{
			
		}
	}

	
}
