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
	public String htmlChart(boolean strandSpecific, boolean normalizedByCounts, boolean range0to1) 
	throws Exception {

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

//		plot.addShapeMarkers(Shape.DIAMOND, Color.BLACK, 4);
		plot.setColor(c);
		
		return plot;
	}
	
	private Data makeData(FeatAligner.PlotType type, boolean normalizedByCounts, boolean range0to1)
	throws Exception
	{
		int NUMPOINTS = 500;

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

		double[] plotArrSmall = plotArr;
		if (plotArrSmall.length > NUMPOINTS)
		{
			plotArrSmall = new double[NUMPOINTS];
			MatUtils.downscaleArray(plotArrSmall, plotArr);
		}
		MatUtils.nansToVal(plotArrSmall, -1.0);
		
//		System.err.println("type= " + type + ", totalsArr=" + ListUtils.excelLine(totalsArr));
//		System.err.println("type= " + type + ", countsArr=" + ListUtils.excelLine(countsArr));
//		System.err.println("type= " + type + ", plotArrSmall=" + ListUtils.excelLine(plotArr));

		
		Data data;
		if (range0to1)
		{
			lastMin = 0.0;
			lastMax = 1.0;
		}
		else
		{
			lastMin = 0.0; // StatUtils.min(plotArrSmall);
			lastMax = StatUtils.max(plotArrSmall);
		}
		data = DataUtil.scaleWithinRange(lastMin, lastMax, plotArrSmall);

		return data;
	}

	@Override
	public void addAlignmentPos(int genomeRelPos, double fwStrandScore,
			double revStrandScore, String featName, String featChr,
			int featCoord, Strand featStrand) {
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
