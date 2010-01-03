/**
 * 
 */
package edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.List;

import javax.swing.text.html.HTML;

import com.googlecode.charts4j.Color;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;


/**
 * @author benb
 *
 * The client streams CpGs into the streamCpg function, and the summarizer can
 * give various statistics.
 * 
 */
public abstract class CpgSummarizer {
	
	protected MethylDbQuerier querier = null;
//	protected String sampleName = "UNKNOWN_SAMPLE";
//	protected String featName = "UNKNOWN_FEAT";
	protected int nCpgsSeen = 0;
	
	// Arbitrary summary statistics for a single value.
	protected double numVals = 0.0;
	protected double valsTotal = 0.0;
	protected double valsSquareTotal = 0.0;
	protected double valsMin = Double.NaN;
	protected double valsMax = Double.NaN;
	
	// Used for chart drawing
	protected String label = "UNKNOWN_SUMMARIZER";
	protected String desc = "NO_DESC";
	protected double naturalMin = Double.NaN;
	protected double naturalMax = Double.NaN;
	protected Color color = Color.BLACK;
	
	/**
	 * We need the querier criteria for certain output
	 */
	public CpgSummarizer() {
		super();
		init(new MethylDbQuerier());
	}

	
	
	/**
	 * We need the querier criteria for certain output
	 */
	public CpgSummarizer(MethylDbQuerier inQuerier) {
		super();
		init(inQuerier);
	}

//	/**
//	 * We need the querier criteria for certain output
//	 */
//	public CpgSummarizer(MethylDbQuerier inQuerier, String sampleName, String featName) {
//		super();
//		this.setFeatName(featName);
//		this.setSampleName(sampleName);
//		init(inQuerier);
//	}

	// This should be called by all constructors
	protected void init(MethylDbQuerier inQuerier)
	{
		this.querier = inQuerier;
		
		this.setDesc(String.format("minReads=%d nonconvFilt=%s, maxOppAfrac=%f",
				querier.getMinCTreads(), ""+querier.getUseNonconversionFilter(), querier.getMaxOppstrandAfrac()));

	}
	
	public void reset()
	{
		nCpgsSeen = 0;
		numVals = 0.0;
		valsTotal = 0.0;
		valsSquareTotal = 0.0;
		valsMin = Double.POSITIVE_INFINITY;
		valsMax = Double.POSITIVE_INFINITY;
	}
	
	
	/**
	 * Override this!
	 * @param cpg
	 */
	public void streamCpg(Cpg cpg)
	{
		nCpgsSeen++;
	}
	
	
		
	
	/**
	 * Use this to update summary stats
	 * @param val
	 */
	protected void streamValue(double val)
	{
		if (!Double.isNaN(val))
		{
			this.numVals++;
			this.valsTotal += val;
			this.valsSquareTotal += Math.pow(val, 2.0);
			if (Double.isNaN(this.valsMin) || (val < this.valsMin)) this.valsMin = val;
			if (Double.isNaN(this.valsMax) || (val > this.valsMax)) this.valsMax = val;
		}
	}

	
	
	public String getDesc() {
		return desc;
	}



	public void setDesc(String desc) {
		this.desc = desc;
	}



	public String getSampleName() {
		return this.querier.getMethylTablePrefix();
	}

//	public void setSampleName(String sampleName) {
//		this.sampleName = sampleName;
//	}

	public String getFeatName() {
		return querier.getFeatFilterTypesList();
	}

//	public void setFeatName(String featName) {
//		this.featName = featName;
//	}

	public int getnCpgsSeen() {
		return nCpgsSeen;
	}

	public double getValMean()
	{
		return this.valsTotal / this.numVals;
	}
	
	public double getValStdev()
	{
		return Math.sqrt((valsSquareTotal/numVals)-Math.pow(this.getValMean(),2));
	}
	
	public double getValsTotal() {
		return valsTotal;
	}

	public double getValsSquareTotal() {
		return valsSquareTotal;
	}

	public double getValsMin() {
		return valsMin;
	}

	public double getValsMax() {
		return valsMax;
	}

	public double getNaturalMin() {
		return naturalMin;
	}

	public double getNaturalMax() {
		return naturalMax;
	}

	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	public Color getColor() {
		return color;
	}

	public void setColor(Color color) {
		this.color = color;
	}

	public void setNaturalMin(double naturalMin) {
		this.naturalMin = naturalMin;
	}

	public void setNaturalMax(double naturalMax) {
		this.naturalMax = naturalMax;
	}
	

	/**** HTML tables , Static ****/
	
	public static String htmlTableStart()
	{
		StringWriter sw = new StringWriter(500);
		PrintWriter pw = new PrintWriter(sw);
		
		pw.print("<TABLE BORDER=1>");
		return sw.toString();
	}
		
	public static String htmlTableFinish()
	{
		StringWriter sw = new StringWriter(500);
		PrintWriter pw = new PrintWriter(sw);
		
		pw.print("</TABLE>");
		return sw.toString();
	}

	public static String htmlTable(List<CpgSummarizer> summarizers)
	{
		StringWriter sw = new StringWriter(100000);
		PrintWriter pw = new PrintWriter(sw);
		
		pw.print(htmlTableStart());
		for (CpgSummarizer summ : summarizers)
		{
			pw.print(summ.htmlTableRow());
		}
		pw.print(htmlTableFinish());
		
		return sw.toString();
	}
	

	/**** HTML tables , non-static ****/
	public String htmlTableRow()
	{
		StringWriter sw = new StringWriter(30000);
		PrintWriter pw = new PrintWriter(sw);

		pw.println("<TR>");
		pw.printf("<TD>%s</TD>\n", this.getLabel());
		pw.printf("<TD>%s</TD>\n", this.getSampleName());
		pw.printf("<TD>%s</TD>\n", this.getFeatName());
		pw.printf("<TD>%d</TD>\n", this.getnCpgsSeen());
		pw.printf("<TD>%.3f</TD>\n", this.getValMean());
		pw.printf("<TD>%.3f</TD>\n", this.getValStdev());
		pw.printf("<TD>%s</TD>\n", this.getDesc());
		pw.print("</TR>");
		
		
		return sw.toString();
	}
	
	

}