package edu.usc.epigenome.genomeLibs.FeatAligners;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

import org.apache.commons.math.stat.StatUtils;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;

import com.googlecode.charts4j.*;

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;

public class FeatAlignerAveraging extends FeatAligner {

	// arr[0] fwTotalScores, arr[1] fwTotalFeats, arr[2] revTotalScores, arr[3] revTotalFeats
	double[][] arr;
	Set<GenomicRange> featsSeen;

	
	protected final int MAX_NUMPOINTS = 500;
	protected final int NUMPOINTS_STEP = 2;
	protected final int MINCOUNTS = 100; // Set to 0 for readCounts
//	double MIN_FRACTION_PASSING = 0.5;
	protected final double MIN_FRACTION_CONTIGUOUS_STRETCH = 0.3;

	
	double lastMin;
	double lastMax;
	
	/**
	 * @param flankSize
	 * @param zeroInit
	 */
	public FeatAlignerAveraging(int flankSize, boolean zeroInit) {
		super(flankSize, zeroInit);

		int nC = (flankSize*2) + 1;
		this.arr = new double[4][nC];
		
		MatUtils.initMat(arr, (zeroInit)  ? 0.0 : Double.NaN);
		
		// We want absolute numbers, so count the number of feats
		featsSeen = new HashSet<GenomicRange>();
	}

	
	
	@Override
	public int numCols() {
		return arr[0].length;
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
			if (Double.isNaN(this.arr[0][relPos]))
			{
				this.arr[0][relPos] = fwScore;
				this.arr[1][relPos] = 1.0;
			}
			else
			{
				this.arr[0][relPos] += fwScore;
				this.arr[1][relPos]++;
			}
		}

		if (!Double.isNaN(revScore))
		{
			if (Double.isNaN(this.arr[2][relPos]))
			{
				this.arr[2][relPos] = revScore;
				this.arr[3][relPos] = 1.0;
			}
			else
			{
				this.arr[2][relPos] += revScore;
				this.arr[3][relPos]++;
			}
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
	public String htmlChart(boolean strandSpecific, boolean normalizedByCounts, boolean range0to1, String sample, String feature) 
	throws Exception {

		
		// sorry, special cases
		System.err.println("Got sample" + sample);
		sample = sample.replaceAll("methylCGsRich_", "");
		sample = sample.replaceAll("_", "");
		feature = feature.replaceAll("wgEncodeBroadChipSeqPeaks", "");
		feature = feature.replaceAll(".broadPeak.hg18.nodups.gtf", "");
		feature = feature.replaceAll("xie.*_Encode", "");
		feature = feature.replaceAll(".hg18.gtf", "");
		
		StringBuilder sb = new StringBuilder(10000);
		
		try
		{

		List<Plot> plots = new ArrayList<Plot>(2);
		
		if (strandSpecific)
		{
			plots.add(this.makePlot(PlotType.FW, normalizedByCounts, range0to1));
			plots.add(this.makePlot(PlotType.REV, normalizedByCounts, range0to1));
		}
		else
		{
			plots.add(this.makePlot(PlotType.COMBINED, normalizedByCounts, range0to1));
		}

        LineChart chart = GCharts.newLineChart(plots);
        chart.setSize(600, 200);
        chart.setDataEncoding(DataEncoding.SIMPLE);
  
        AxisLabels xAxis = AxisLabelsFactory.newNumericRangeAxisLabels(-this.flankSize, this.flankSize);
        chart.addXAxisLabels(xAxis);
        
        AxisLabels yAxis;
        if (range0to1)
        {
        	yAxis = AxisLabelsFactory.newNumericRangeAxisLabels(0.0, 1.0);
        }
        else
        {
        	double min = this.lastMin;
        	double max = this.lastMax;
        	yAxis = AxisLabelsFactory.newNumericRangeAxisLabels(min, max);
        	System.err.printf("min=%f, max=%f\n", min, max);
        }
        chart.addYAxisLabels(yAxis);
		chart.addRightAxisLabels(yAxis);
       
        chart.setTitle(String.format("feat=%s, me=%s", feature, sample));
        
        sb.append(String.format("<P>strandSpec=<EM>%s</EM> normalizedByCounts=<EM>%s</EM>, range0to1=<EM>%s</EM></P>",
				strandSpecific, normalizedByCounts, range0to1));
        
        sb.append("<IMG ");
        sb.append(" SRC=\"");
        sb.append(chart.toURLString());
        sb.append("\" ALT=\"google chart\">\n");
		}
		catch (Exception e)
		{
			System.err.println("Could not make chart: " + e.toString());
			e.printStackTrace(System.err);
		}

        
        return sb.toString();
	}
	
	
	/**
	 * @param type 1=fw, 2=rev, 3=combined
	 * @param normalizedByCounts
	 * @return
	 */
	private Plot makePlot(FeatAligner.PlotType type, boolean normalizedByCounts, boolean range0to1)
	throws Exception
	{
		Data data = this.makeData(type, normalizedByCounts, range0to1);
		Plot plot = Plots.newPlot(data);
		
		Color c;
		switch (type)
		{
		case FW:
			c = Color.BLUE; break;
		case REV:
			c = Color.RED; break;
		default:
			c = Color.CORNFLOWERBLUE; break;
		}

		// Bug in charts4j makes it so this doesn't work with 
		// charts with lots of data points.  There is a proposed 
		// workaround which i added to the project.
		if (data.getSize() < 250)
		{
			plot.addShapeMarkers(Shape.DIAMOND, Color.BLACK, 4);
		}
		
		// Add a line halfway
		int midpoint = (int)Math.round((double)data.getSize()/2.0);
		plot.addShapeMarker(Shape.VERTICAL_LINE_FULL, Color.BLACK,1,midpoint);
		
		plot.setColor(c);
		
		return plot;
	}
	
	private Data makeData(FeatAligner.PlotType type, boolean normalizedByCounts, boolean range0to1)
	throws Exception
	{

		double[] totalsArr;
		double[] countsArr;
		switch (type)
		{
		case FW:
			totalsArr = this.arr[0];
			countsArr = this.arr[1];
			break;
		case REV:
			totalsArr = this.arr[2];
			countsArr = this.arr[3];
			break;
		case COMBINED:
			totalsArr = MatUtils.vectSum(this.arr[0], this.arr[2]);
			countsArr = MatUtils.vectSum(this.arr[1], this.arr[3]);
			break;
		default:
			totalsArr = null; countsArr = null;
			break;
		}
		
		// Normalize
		double[] plotArr = totalsArr;
		if (normalizedByCounts)
		{
			plotArr = MatUtils.divVects(totalsArr, countsArr);
		}
		else
		{
			// Always divide by num feats
			plotArr = MatUtils.divVect(totalsArr, (double)this.numFeats());
		}

		// Do lower resolution until we get enough bins with high enough numbers
		// of counts.  Check 
		int numpointsCur = MAX_NUMPOINTS;
		double fracPassing = 0.0;
		double maxFracContiguous = 0.0;
		double[] plotArrSmall = null;
		while ((numpointsCur > 0) && (maxFracContiguous < MIN_FRACTION_CONTIGUOUS_STRETCH))
		{
			plotArrSmall = plotArr;
			if (plotArrSmall.length > numpointsCur)
			{
				plotArrSmall = new double[numpointsCur];
				MatUtils.downscaleArray(plotArrSmall, plotArr);
			}

			// Impose a minimum number of points to display
			double[] avgCounts = new double[numpointsCur];
			MatUtils.downscaleArray(avgCounts, countsArr);
			fracPassing = 0.0;
			maxFracContiguous = 0.0;
			int beginStretch = 0;
			for(int i = 0; i < numpointsCur; i++)
			{
				double totalCounts = (double)avgCounts[i] * ((double)countsArr.length/(double)numpointsCur);
				if (totalCounts >= MINCOUNTS)
				{
					fracPassing += 1.0/(double)numpointsCur;
					
					int numContig = i - beginStretch;
					maxFracContiguous = Math.max(maxFracContiguous, (double)numContig/(double)numpointsCur);
				}
				else
				{
					beginStretch = i;
					plotArrSmall[i] = Double.NaN;
				}
			}
			
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe(String.format(
					"Checked chart with %d points, found %f passing (maxcontig=%f)\n", numpointsCur, fracPassing, maxFracContiguous));

			// Increment
			numpointsCur -= NUMPOINTS_STEP; 
		}

		
		
//		System.err.println("type= " + type + ", totalsArr=" + ListUtils.excelLine(totalsArr));
//		System.err.println("type= " + type + ", countsArr=" + ListUtils.excelLine(countsArr));
//		System.err.println("type= " + type + ", plotArr=" + ListUtils.excelLine(plotArr));
//		System.err.println("type= " + type + ", plotArr=" + ListUtils.excelLine(plotArrSmall));

		
		if (range0to1)
		{
			lastMin = 0.0;
			lastMax = 1.0;
		}
		else
		{
			lastMin = ((MatUtils.nanMin(plotArrSmall)) > 0.0) ? 0.0 : MatUtils.nanMin(plotArrSmall);
			lastMax = MatUtils.nanMax(plotArrSmall);
			System.err.println("min="+lastMin+"\tmax="+lastMax);
		}

		// If it's still nan, we had no values!
		if (Double.isNaN(lastMin))
		{
			lastMin = 0.0;
			lastMax = 1.0;
		}
		
		
		// Not allowed to have NANs in array
		MatUtils.nansToVal(plotArrSmall, lastMin-1.0);
		Data data = DataUtil.scaleWithinRange(lastMin, lastMax, plotArrSmall);

		return data;
	}

	@Override
	public void addAlignmentPos(int genomeRelPos, double fwStrandScore,
			double revStrandScore, String featName, String featChr,
			int featCoord, Strand featStrand, double sortVal) {
		int colInd = this.getColumnInd(genomeRelPos, featCoord, featStrand);
		
		this.addAlignmentPosFast(colInd,
				(featStrand == StrandedFeature.NEGATIVE) ? revStrandScore : fwStrandScore,
						(featStrand == StrandedFeature.NEGATIVE) ? fwStrandScore : revStrandScore);

		GenomicRange gr = new GenomicRange(featChr, featCoord, featCoord, featStrand);
		featsSeen.add(gr);

	}

	public int numFeats()
	{
		return featsSeen.size();
	}
	
	public void matlabCsv(PrintWriter pw, boolean strandSpecific)
	{
		if (strandSpecific)
		{
			if (this.zeroInit)
			{
				// If it's zero init, it means we just want to divide by the
				// total.
				MatUtils.initMat(this.arr[1],(double)this.numFeats());
				MatUtils.initMat(this.arr[3],(double)this.numFeats());
			}
			MatUtils.matlabCsv(pw, this.arr);
		}
		else
		{
			
		}
	}

	
}
