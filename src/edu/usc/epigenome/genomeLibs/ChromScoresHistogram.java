package edu.usc.epigenome.genomeLibs;

import java.awt.Dimension;
import java.awt.GradientPaint;

import java.awt.Color;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.xy.*;
import org.jfree.data.xy.*;
import org.jfree.ui.ApplicationFrame;
import org.jfree.data.statistics.*;

public class ChromScoresHistogram extends ApplicationFrame {

	static final long serialVersionUID = 1000111;
	
	public ChromScoresHistogram(String title, double []series, int num_bars) {
		
		super(title);
		
		System.err.println("Creating dataset");
		
		IntervalXYDataset dataset = createDataset(series,num_bars);
		
		JFreeChart chart = createHistChart(dataset, title);
		
		// Make the panel
		ChartPanel chartPanel = new ChartPanel(chart, false);
		chartPanel.setPreferredSize(new Dimension(500, 270));
		setContentPane(chartPanel);
	
	}



	private static JFreeChart createHistChart(IntervalXYDataset dataset, String title) {

//		create the chart...
		JFreeChart chart = ChartFactory.createHistogram(
				title,       // chart title
				"Category",               // domain axis label
				"Value",                  // range axis label
				dataset,                  // data
				PlotOrientation.VERTICAL, // orientation
				true,                     // include legend
				true,                     // tooltips?
				false                     // URLs?
		);


//		NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...

//		set the background color for the chart...
		chart.setBackgroundPaint(Color.white);

//		get a reference to the plot for further customisation...
		XYPlot plot = (XYPlot) chart.getPlot();
		plot.setBackgroundPaint(Color.lightGray);
		plot.setDomainGridlinePaint(Color.white);
		plot.setDomainGridlinesVisible(true);
		plot.setRangeGridlinePaint(Color.white);


//		set the range axis to display integers only...
		final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
		rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());

//		disable bar outlines...
		XYBarRenderer renderer = (XYBarRenderer) plot.getRenderer();
		
		renderer.setMargin(0.2);

//		set up gradient paints for series...
		GradientPaint gp0 = new GradientPaint(0.0f, 0.0f, Color.blue, 
				0.0f, 0.0f, new Color(0, 0, 64));
		GradientPaint gp1 = new GradientPaint(0.0f, 0.0f, Color.green, 
				0.0f, 0.0f, new Color(0, 64, 0));
		GradientPaint gp2 = new GradientPaint(0.0f, 0.0f, Color.red, 
				0.0f, 0.0f, new Color(64, 0, 0));
		renderer.setSeriesPaint(0, gp0);
		renderer.setSeriesPaint(1, gp1);
		renderer.setSeriesPaint(2, gp2);

				
//		CategoryAxis domainAxis = plot.getDomainAxis();
//		domainAxis.setCategoryLabelPositions(
//				CategoryLabelPositions.createUpRotationLabelPositions(
//						Math.PI / 6.0));
////		OPTIONAL CUSTOMISATION COMPLETED.

		return chart;

	}


	private static IntervalXYDataset createDataset(double[] series, int num_bars) {

		HistogramDataset dataset = new HistogramDataset();
		dataset.addSeries(new Double(1), series, num_bars);
		return dataset;
	}

	//  DEMO DEMO
//	private static IntervalXYDataset createDataset(ChromScoresFast csf, int num_bars) {
//
//		Random rand = new Random();
//		
//		int n = 100000;
//		double[] series1 = new double[n];
//		double[] series2 = new double[n];
//		for (int i = 0; i < n; i++)
//		{
//			series1[i] = rand.nextGaussian();
//			series2[i] = 5 + rand.nextGaussian();
//		}
//		
//		
//		// create the dataset...
//		HistogramDataset dataset = new HistogramDataset();
//		dataset.addSeries(new Double(1), series1, 100);
//		dataset.addSeries(new Double(1), series2, 100);
//
//
//		return dataset;
//	}	
	
	
}
