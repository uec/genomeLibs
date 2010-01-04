package edu.usc.epigenome.genomeLibs.FeatAligners;

import heatMap.Gradient;
import heatMap.HeatMap;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
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

	public List<Color> ColorCycle = Arrays.asList(Color.DARKBLUE, Color.BLUE, Color.BLUEVIOLET, Color.VIOLET, Color.LIGHTSALMON, Color.RED);
	
	// i = type: arr[0] fwTotalScores, arr[1] revTotalScores, 
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
		int nEls = 2*nFeats*nC;
		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("About to allocate memory ("+nEls+" doubles) for alignerEachfeat()\n");
		
		if ((nEls*4) > 1E9)
		{
			System.err.printf("Trying to allocated an array of %d elements (2 strands, %d cols, %d feats) - too big\n",
					nEls, nC, nFeats);
			System.exit(1);
		}
		
		// i = type: arr[0] fwTotalScores, arr[1] revTotalScores
		// j = featNum: arr[0][5] = fwTotalScores for feat 6
		// k = coordinate: arr[0][5][350] = coord (350 - flank) relative to feature center.
		this.arr = new double[2][nFeats][nC];
		this.featCoords = new GenomicRange[nFeats];
		this.featNames = new String[nFeats];
		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("About to initMat()\n");
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

			double[][] data = MatUtils.nanMeanMats(this.arr[0], this.arr[1]);
			
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(String.format(
					"num nans: this.arr[0]=%d, this.arr[1]=%d, sum(this.arr[0..1])=%d\n",
					MatUtils.countNans(this.arr[0]), MatUtils.countNans(this.arr[1]), MatUtils.countNans(data)));
		
//			data = new double[10][30];
//			MatUtils.initMat(data, Double.NaN);
//			data[0][5] = 5;
//			data[3][0] = 3;
			

			data = MatUtils.sortRows(data,-0.3333,10);
			data = MatUtils.downscaleMatRows(data, 200, 4.0);
			data = MatUtils.downscaleMatCols(data,5, 0.0);

			// NO NANs allowed
			MatUtils.nansToVal(data, 0.0);

			List<Data> chartData = DataUtil.scale(data);
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
		if (strandSpecific)
		{
			MatUtils.matlabCsv(pw, this.arr[0], this.numFeats(), 0);
			MatUtils.matlabCsv(pw, this.arr[1], this.numFeats(), 0);
		}
		else
		{
			
		}
	}
	
	public HeatMap heatMap(double[] colorMinMax)
	throws Exception
	{
		double[][] data = MatUtils.nanMeanMats(this.arr[0], this.arr[1]);
		
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(String.format(
				"num nans: this.arr[0]=%d, this.arr[1]=%d, sum(this.arr[0..1])=%d\n",
				MatUtils.countNans(this.arr[0]), MatUtils.countNans(this.arr[1]), MatUtils.countNans(data)));
	
//		data = new double[10][30];
//		MatUtils.initMat(data, Double.NaN);
//		data[0][5] = 5;
//		data[3][0] = 3;
		

		data = MatUtils.sortRows(data,-0.3333,10);
		data = MatUtils.downscaleMatRows(data, 500, 20.0);
		data = MatUtils.downscaleMatCols(data, 30, 2.0);

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
