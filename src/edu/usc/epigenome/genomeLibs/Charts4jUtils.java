package edu.usc.epigenome.genomeLibs;

import edu.usc.epigenome.genomeLibs.MatUtils;

import com.googlecode.charts4j.AxisLabels;
import com.googlecode.charts4j.AxisLabelsFactory;
import com.googlecode.charts4j.BarChart;
import com.googlecode.charts4j.BarChartPlot;
import com.googlecode.charts4j.Color;
import com.googlecode.charts4j.Data;
import com.googlecode.charts4j.DataEncoding;
import com.googlecode.charts4j.DataUtil;
import com.googlecode.charts4j.GCharts;
import com.googlecode.charts4j.Plots;
import com.googlecode.charts4j.Shape;

public class Charts4jUtils {

	
	public static String histogramBarPlot(double[] vals, double min, double max, int nBins, String name, Color color)
	{
		// First we get the bin counts
		double[] binVals = MatUtils.histogramBins(vals, min, max, nBins, true);
		double valMean = MatUtils.nanMean(vals);
		
		final double Y_MIN = 0;
		final double Y_MAX = 0.5;
		
		// Nan breaks charts4j
		MatUtils.nansToVal(binVals, -0.1);
		
		//Clip those over YMAX, it doesn't deal well with them.
		MatUtils.clipToMax(binVals, Y_MAX);
		
		Data data = DataUtil.scaleWithinRange(Y_MIN, Y_MAX, binVals); // 0-1 because it's normalized (maybe it
		BarChartPlot plot = Plots.newBarChartPlot(data, color);
	
		// Make a vertical line at the mean
		double meanRel = (valMean-min)/(max-min);
		int linePos = (int)Math.round((double)nBins*meanRel);
		try
		{
			plot.addShapeMarker(Shape.VERTICAL_LINE_FULL, Color.BLACK,1,linePos);
		}
		catch (Exception e)
		{
			System.err.println("Couldn't draw mean line:\n" + e.toString());
			e.printStackTrace();
		}
	
		
		BarChart chart = GCharts.newBarChart(plot);
		chart.setSize(200, 130); // Make sure to make it wide enough, it doesn't autoscale
        chart.setDataEncoding(DataEncoding.EXTENDED);

        AxisLabels xAxis = AxisLabelsFactory.newNumericRangeAxisLabels(min, max, 2*((min-max)/(double)nBins));
       	chart.addXAxisLabels(xAxis);
//       	chart.addXAxisLabels(AxisLabelsFactory.newAxisLabels("Features"));
       	

       	
        //AxisLabels yAxis = AxisLabelsFactory.newNumericRangeAxisLabels(Y_MIN, Y_MAX); // Make them %
        AxisLabels yAxis = AxisLabelsFactory.newNumericRangeAxisLabels(Y_MIN*(double)100, Y_MAX*(double)100, 20.0);
       	chart.addYAxisLabels(yAxis);
       	chart.addRightAxisLabels(yAxis);

       	chart.setTitle(name);
       	chart.setSpaceBetweenGroupsOfBars(1);
       	chart.setSpaceWithinGroupsOfBars(1);
       	chart.setBarWidth(10);
       	chart.setLegendMargins(0, 0);
       	
       	
        StringBuffer sb = new StringBuffer(10000);
 		sb.append("<IMG ");
		sb.append(" SRC=\"");
		sb.append(chart.toURLString());
		sb.append("\" ALT=\"google chart\">\n");
		
		return sb.toString();
	}
	
}
