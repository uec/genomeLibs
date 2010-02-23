package edu.usc.epigenome.genomeLibs.FeatAligners;

import heatMap.Gradient;
import heatMap.HeatMap;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Logger;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.seq.StrandedFeature.Strand;

import com.googlecode.charts4j.AxisLabels;
import com.googlecode.charts4j.AxisLabelsFactory;
import com.googlecode.charts4j.Color;
import com.googlecode.charts4j.Data;
import com.googlecode.charts4j.DataEncoding;
import com.googlecode.charts4j.DataUtil;
import com.googlecode.charts4j.GCharts;
import com.googlecode.charts4j.LineChart;
import com.googlecode.charts4j.Plot;
import com.googlecode.charts4j.Plots;

import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAligner.PlotType;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;

public class FeatAlignerEachfeat extends FeatAligner {

	public List<Color> ColorCycle = Arrays.asList(Color.DARKBLUE, Color.BLUE, Color.BLUEVIOLET, Color.VIOLET, Color.PINK, Color.LIGHTSALMON, Color.RED);
	
	protected final static int NBINS = 2;
	protected final static int BP_SMOOTHING = 400;
	protected final static int HEATMAP_ROWS = 100;
	
	// i = type: arr[0] fwTotalScores, arr[1] revTotalScores, arr[2] fwNumScores, arr[3] revNumScores,
	// j = featNum: arr[0][5] = fwTotalScores for feat 6
	// k = coordinate: arr[0][5][350] = coord (350 - flank) relative to feature center.
	protected double[][][] arr;
	protected String[] featNames;
	protected Double[] sortVals;
	protected GenomicRange[] featCoords;
	protected Map<GenomicRange,Integer> featinds = new TreeMap<GenomicRange,Integer>();
	protected int nFeatsSeen = 0;
	
	/**
	 * @param flankSize
	 * @param zeroInit
	 * @param nFeats: If it's more than the actual number, it's ok
	 * @param downsampleCols: If it's zero, we just use (2*flankSize)+1
	 */
	public FeatAlignerEachfeat(int inFlankSize, boolean zeroInit, int nFeats, int inDownscaleCols) {
		super(inFlankSize, zeroInit);

		this.flankSize = inFlankSize;
		this.downscaleCols = (inDownscaleCols>0) ? inDownscaleCols : ((flankSize*2)+1); 
		int nEls = 2*nFeats*this.downscaleCols;
		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("About to allocate memory ("+nEls+" doubles) for alignerEachfeat() (" + nFeats + ") features\n");
		
		if ((nEls*4) > 1E9)
		{
			System.err.printf("Trying to allocated an array of %d elements (2 strands, %d cols, %d feats) - too big\n",
					nEls, this.downscaleCols, nFeats);
			System.exit(1);
		}
		
		// i = type: arr[0] fwTotalScores, arr[1] revTotalScores, arr[2] fwNumScores, arr[3] revNumScores,
		// j = featNum: arr[0][5] = fwTotalScores for feat 6
		// k = coordinate: arr[0][5][350] = coord (350 - flank) relative to feature center.
		this.arr = new double[4][nFeats][this.downscaleCols];
		this.featCoords = new GenomicRange[nFeats];
		this.featNames = new String[nFeats];
		this.sortVals = new Double[nFeats];
		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("About to initMat()\n");
		MatUtils.initMat(arr, (zeroInit)  ? 0.0 : Double.NaN);
	}

	@Override
	public int numCols() {
		return arr[0][0].length;
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
			String featName, String featChr, int featCoord, Strand featStrand, double sortVal) {

		GenomicRange gr = new GenomicRange(featChr, featCoord, featCoord, featStrand);
		int featInd = this.getInd(gr, featName);
		int colInd = this.getColumnInd(genomeRelPos, featCoord, featStrand, true);
		this.sortVals[featInd] = new Double(sortVal);
		
		// Flip the strand of the scores if features are flipped
		double fwScoreRel = (featStrand==StrandedFeature.NEGATIVE) ? revStrandScore : fwStrandScore;
		double revScoreRel = (featStrand==StrandedFeature.NEGATIVE) ? fwStrandScore : revStrandScore;

		// Add it to old score
		arr[0][featInd][colInd] = MatUtils.nanSum(fwScoreRel, arr[0][featInd][colInd]);
		arr[1][featInd][colInd] = MatUtils.nanSum(revScoreRel, arr[1][featInd][colInd]);

		// And add to totals
		if (!Double.isNaN(fwScoreRel)) arr[2][featInd][colInd] = MatUtils.nanSum(1.0, arr[2][featInd][colInd]);
		if (!Double.isNaN(revScoreRel))  arr[3][featInd][colInd] = MatUtils.nanSum(1.0, arr[3][featInd][colInd]);
		
	}



	public Double[] sortRowsExponential(double exponentialFactor, int minDist)
	{
		int nCols = this.arr[0][0].length;
		return this.sortRowsExponential(exponentialFactor, minDist, 0, nCols-1);
	}
	
	public Double[] sortRowsExponential(double exponentialFactor, int minDist, int colStart, int colEnd)
	{
		double[][] dataFull = MatUtils.nanMeanMats(this.arr[0], this.arr[1]);
		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(
				String.format("Sorting %d rows EXPONENTIAL %f, %d\n",this.nFeatsSeen, exponentialFactor, minDist));
		// Did we actually see as many features as expected?
		Double[] sortVals = new Double[this.nFeatsSeen];
		for (int i = 0; i < this.nFeatsSeen; i++)
		{
			double mean = (exponentialFactor!=0.0) ? 
					MatUtils.nanMeanExponentialWeighting(dataFull[i], exponentialFactor,minDist,colStart,colEnd) : 
						MatUtils.nanMean(dataFull[i], colStart, colEnd);	
					sortVals[i] = new Double(mean);
		}
			
		
		this.sortRowsByList(sortVals);
		return sortVals;
	}

	public void sortRowsByList(Double[] sortVals)
	{
		for (int i = 0; i < this.arr.length; i++)
		{
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(
					String.format("Sorting %d rows by list arr[%d]\n",this.nFeatsSeen,i));
			this.arr[i] = MatUtils.sortRowsByList(this.arr[i], sortVals);
		}
	}
	
	public Double[] sortRowsBySortVals()
	{
		// Only those that we've actually seen.
		Double[] sortList = Arrays.copyOf(this.sortVals, this.nFeatsSeen);
		
		sortRowsByList(sortList);
		return sortList;
	}
	
	
//	@Override
//	public String htmlChart(boolean strandSpecific, boolean normalizedByCounts, boolean range0to1) throws Exception{
//		FeatAlignerAveraging av = this.toAverageFeatAligner();
//		return av.htmlChart(strandSpecific, normalizedByCounts, range0to1);

	@Override
	public String htmlChart(boolean strandSpecific, boolean normalizedByCounts, boolean range0to1, String sample, String feature) throws Exception
	{

		StringBuilder sb = new StringBuilder(10000);

		try
		{

			// DOWNSCALING
		
			double[][] fw = MatUtils.divMats(this.arr[0], this.arr[2]);
			double[][] rev = MatUtils.divMats(this.arr[1], this.arr[3]);
			
			double[][] dataFull = MatUtils.nanMeanMats(fw, rev);
			double[][]data = new double[this.nFeatsSeen][];
			for (int i = 0; i < this.nFeatsSeen; i++) data[i] = dataFull[i];
			
			
			// Do rows first since cols has to transpose 
			int downsampleTo = Math.min(this.downscaleCols, 1000); // Harder to go above 500
			double smoothFact = Math.ceil(((2.0*(double)this.flankSize)/(double)this.downscaleCols)/(2.0*(double)BP_SMOOTHING));
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(String.format("Downscaling %d cols to %d (smoothing %f)\n",
					data[0].length,downsampleTo, smoothFact));
			smoothFact = 0.0;
//			data = MatUtils.downscaleMatRows(data, downsampleTo, smoothFact);
			data = MatUtils.downscaleMatRows(data, downsampleTo, 2);
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(String.format("Downscaling %d rows to 5\n",data.length));
			data = MatUtils.downscaleMatCols(data,NBINS, 0.0);

			// NO NANs allowed
			MatUtils.nansToVal(data, 0.0);

			List<Data> chartData = new ArrayList<Data>(data.length);
			for (int i = 0; i < data.length; i++)
			{
				Data newData = (range0to1) ? DataUtil.scaleWithinRange(0.0, 1.0, data[i]) : DataUtil.scale(data[i]);
				chartData.add(newData);
			}
			List<Plot> plots = new ArrayList<Plot>(chartData.size());
			for (int i = 0; i < chartData.size(); i++)
			{
				Data d = chartData.get(i);
				Plot p = Plots.newPlot(d);
				p.setColor(ColorCycle.get( i%ColorCycle.size()));
				plots.add(p);
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
				double min = 0;
				double max = 1;
				yAxis = AxisLabelsFactory.newNumericRangeAxisLabels(min, max);
				System.err.printf("min=%f, max=%f\n", min, max);
			}
			chart.addYAxisLabels(yAxis);
			chart.addRightAxisLabels(yAxis);

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


		return sb.toString();	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.FeatAligners.FeatAligner#toAverageFeatAligner()
	 */
	@Override
	public FeatAlignerAveraging toAverageFeatAligner() 
	{

		FeatAlignerAveraging out = new FeatAlignerAveraging(this.flankSize, this.zeroInit);
		
		for (int j = 0; j < this.numFeats(); j++)
		{
			double[] fwScores = arr[0][j];
			double[] revScores = arr[1][j];
			try
			{
			fwScores = MatUtils.divVects(arr[0][j], arr[2][j]);
			revScores = MatUtils.divVects(arr[1][j], arr[3][j]);
			}
			catch (Exception e)
			{
				System.err.println(e);
				e.printStackTrace();
			}
			out.addAlignmentFast(fwScores, revScores);
		}
		
		return out;
	}

	
	private int getInd(GenomicRange gr, String featName) {
		
		Integer ind = featinds.get(gr);
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(
				String.format("Adding new feature at pos=%d to ind=%d\n",gr.getStart(),ind));
		if (ind == null)
		{
			ind = new Integer(nFeatsSeen++);
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(
					String.format("Adding new feature at pos=%d to ind=%d\n",gr.getStart(),ind));
			featinds.put(gr, ind);
			this.featNames[ind.intValue()] = featName;
			this.featCoords[ind] = gr;
		}
		return ind;
	}

	public void matlabCsv(PrintWriter pw, boolean strandSpecific)
	{
		double[][] fw = this.arr[0];
		double[][] rev = this.arr[1];
		if (!this.zeroInit)
		{
			try 
			{
				fw = MatUtils.divMats(this.arr[0], this.arr[2]);
				rev = MatUtils.divMats(this.arr[1], this.arr[3]);
			}
			catch (Exception e)
			{
				System.err.println(e);
				e.printStackTrace();
			}
		}
		
			
		if (strandSpecific)
		{
			MatUtils.matlabCsv(pw, fw, this.numFeats(), 0);
			MatUtils.matlabCsv(pw, rev, this.numFeats(), 0);
		}
		else
		{
			double[][] dataFull = MatUtils.nanMeanMats(fw, rev);
			System.err.println("Making matlab for " + dataFull.length + " rows (numFeats=" + this.numFeats());
			MatUtils.matlabCsv(pw, dataFull, this.numFeats(), 0);
		}
	}
	
	public HeatMap heatMap(double[] colorMinMax)
	throws Exception
	{
		double[][] fw = this.arr[0];
		double[][] rev = this.arr[1];
		try 
		{
			fw = MatUtils.divMats(this.arr[0], this.arr[2]);
			rev = MatUtils.divMats(this.arr[1], this.arr[3]);
		}
		catch (Exception e)
		{
			System.err.println(e);
			e.printStackTrace();
		}

		double[][] dataFull = MatUtils.nanMeanMats(fw, rev);
		
		// Did we actually see as many features as expected?
		double[][]data = new double[this.nFeatsSeen][];
		for (int i = 0; i < this.nFeatsSeen; i++) data[i] = dataFull[i];

		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(String.format(
				"num nans: this.arr[0]=%d, this.arr[1]=%d, sum(this.arr[0..1])=%d\n",
				MatUtils.countNans(this.arr[0]), MatUtils.countNans(this.arr[1]), MatUtils.countNans(data)));
	
//		data = new double[10][30];
//		MatUtils.initMat(data, Double.NaN);
//		data[0][5] = 5;
//		data[3][0] = 3;
		

	//	data = MatUtils.sortRows(data,-0.3333,10);
		int downsampleTo = Math.min(this.downscaleCols, 500); // Harder to go above 500
		double smoothFact = Math.ceil(((2.0*(double)this.flankSize)/(double)this.downscaleCols)/(2.0*(double)BP_SMOOTHING));
//		data = MatUtils.downscaleMatRows(data, downsampleTo,smoothFact);; 
		data = MatUtils.downscaleMatCols(data, HEATMAP_ROWS, 2.0);

		// NO NANs allowed
		MatUtils.nansToVal(data, 0.0);

		// HeatMap always autoscales, so we have to put a stupid fake upper and lower
		// bound in there if we want colors to scale correctly.
		if (colorMinMax!=null)
		{
			data[0][0] = colorMinMax[0];
			data[0][1] = colorMinMax[1];
		}
		
		HeatMap panel = new HeatMap(MatUtils.transposedMat(data), true, Gradient.GRADIENT_HEAT);

		
		// set miscellaneous settings

		panel.setDrawLegend(true);

//		panel.setTitle("features");
//		panel.setDrawTitle(true);

		panel.setXAxisTitle("distance");
		panel.setDrawXAxisTitle(true);

		panel.setYAxisTitle("feature");
		panel.setDrawYAxisTitle(true);

		panel.setCoordinateBounds(-this.flankSize, this.flankSize, 1, this.nFeatsSeen);

		panel.setDrawXTicks(true);
		panel.setDrawYTicks(true);
		
		return panel;
	}
	
	public void launchSwingHeatmap(double[] colorMinMax)
	throws Exception
	{
		final HeatMap hm = this.heatMap(colorMinMax);
		
		SwingUtilities.invokeLater(new Runnable()
		{
			public void run()
			{
				try
				{
					JFrame hmf = new JFrame();
					hmf.getContentPane().add(hm);
					hmf.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
					
					hmf.setSize(500,500);
					hmf.setVisible(true);
				}
				catch (Exception e)
				{
					System.err.println(e);
					e.printStackTrace();
				}
			}
		}
		);
	}
	

	
}
