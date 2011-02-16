package edu.usc.epigenome.scripts;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import sun.tools.tree.ThisExpression;

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
import com.googlecode.charts4j.Shape;

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.PicardUtils;


public class SamFiveprimeConversion_DEPRECATED_DONT_USE {

	static final String USAGE = "SamFiveprimeConversion -relativeFreqs -outputFilePrefix FiveprimeConversion -minConv 1 input1.sam input2.sam ...";

	static final int BUFFERLEN = 1000;
	static StringBuilder sb = new StringBuilder(BUFFERLEN);
	
	static Pattern pairPat = Pattern.compile("([0-9]*)([ACTGNactgn]?)");

	
	/**
	 * @param args
	 */

	@Option(name="-minConv",usage="minimum number of converted cytosines required")
	protected int minConv = 1;
	@Option(name="-numRecs",usage="If set, this does the first n recs (default off)")
	protected int numRecs = -1;
	@Option(name="-outputFilePrefix",usage="Prefix for output files")
	protected String outputFilePrefix = "FiveprimeConv";
	@Option(name="-numCycles",usage="Number of cycles to track")
	protected int numCycles = 100;
	@Option(name="-minMapq",usage="Minimim mapping quality (default 30)")
	protected int minMapq = 30;
//	@Option(name="-outputReads",usage=" Outputs one line per read (default false)")
//	protected boolean outputReads = false;
	@Option(name="-useCpgsToFilter",usage=" Use CpGs and CpHs to filter if true, otherwise just CpHs (default false)")
	protected boolean useCpgsToFilter = false;
	@Option(name="-noFilters",usage=" Don't display filters (default false)")
	protected boolean noFilters = false;
	@Option(name="-relativeFreqs",usage=" If set, we don't do filtering and we instead plot the position relative to all other positions")
	protected boolean relativeFreqs = false;

	// receives other command line parameters than options
	@Argument
	private List<String> stringArgs = new ArrayList<String>();

	
	public static void clearStrBuf()
	{
		sb.delete(0, BUFFERLEN);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception
	{
		new SamFiveprimeConversion_DEPRECATED_DONT_USE().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception {

		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);
			if (stringArgs.size() < 1) throw new CmdLineException("No input file specified");
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			System.err.println(USAGE);
			// print the list of available options
			parser.printUsage(System.err);
			System.err.println();
			return;
		}
		
				
//		if (outputReads) System.out.println(headerLine());
		

		// Setup cycle counter list
		Map<String,ByCycleCounter> cycleCounters = new HashMap<String,ByCycleCounter>();
		
		for (String fn : stringArgs)
		{

			File inputSamOrBamFile = new File(fn);


			final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
			inputSam.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
			int recCounter = 0;
			CloseableIterator<SAMRecord> samIt = inputSam.iterator();
			SAMREC: while ( !((this.numRecs>0) && (recCounter>this.numRecs)) && samIt.hasNext())
			{
			//record: for (final SAMRecord samRecord : inputSam) {
				// Convert read name to upper case.
				//samRecord.setReadName(samRecord.getReadName().toUpperCase());
				SAMRecord samRecord = samIt.next();
				
				int mapQual = samRecord.getMappingQuality();
				boolean unmapped = samRecord.getReadUnmappedFlag();
				if (unmapped || (mapQual < this.minMapq))
				{
					continue SAMREC;
				}

				String seq = PicardUtils.getReadString(samRecord, true);

				recCounter++;
				if ((recCounter % 1E5)==0) 
					System.err.println("On new record #" + recCounter); 

				try
				{

					String ref = PicardUtils.refStr(samRecord, true);
					//System.err.println("\tref=" + ref);

					if (seq.length() != ref.length())
					{
						System.err.println("SeqLen(" + seq.length() + ") != RefLen(" + ref.length() + ")");
						System.err.println(seq + "\n" + ref);
					}

					ByCycleCounter indReadCounter = new ByCycleCounter(numCycles);
					
					int numConverted = 0;
					int convStart = Integer.MAX_VALUE;
					int seqLen = Math.min(seq.length(), ref.length());
					for (int i = 0; i < this.numCycles; i++)
					{
						if (isCytosine(i,ref,seq))
						{
							boolean cpg = isCpg(i,ref,seq);
							boolean conv = isConverted(i,ref,seq);

							if (conv && (useCpgsToFilter || !cpg)) numConverted++;

							if ((convStart==Integer.MAX_VALUE) && (numConverted>=minConv) )
							{
								convStart = i;
							}

							// The first one is ok to use for actual data, but not
							// ok to use for these statistics, because it throws
							// off the ratios (100% of CpGs at the first base with
							// filtering will be converted, which is not useful information)
//							if (this.relativeFreqs || (convStart != i))
//							{
								// Get counter
								String key = cytosineContext(i, ref, seq);
								ByCycleCounter counter = cycleCounters.get(key);
								if (counter == null)
								{
									counter = new ByCycleCounter(numCycles);
									cycleCounters.put(key, counter);
								}

								boolean filter = (i<convStart);
								counter.increment(filter, conv, i);
//							}

							//if (cpg && (i<convStart)) System.err.printf("Rec %d\tpos=%d\n",recCounter,i);

						}
					}
				}
				catch (Exception e)
				{
					System.err.println("-----------------------------------------");
					System.err.println("Couldn't handle seq #" + recCounter);
					System.err.println(seq);
					e.printStackTrace(System.err);
					System.err.println("-----------------------------------------");
				}

			} // record
			
			System.err.printf("File %s found %d reads\n", fn, recCounter);


			inputSam.close();

		} // File
		
		// HTML output
		String htmlFn = String.format("%s.htm", this.outputFilePrefix);
		PrintWriter pw = new PrintWriter(new File(htmlFn));
		pw.printf("<H1>5prime conversion for %d files</H1><P>\n", stringArgs.size());
		for (String fn : stringArgs)
		{
			pw.printf("%s<BR>\n",fn);
		}
		
		if (!this.relativeFreqs)
		{
			pw.print("</P>\n");
			pw.print("<H4>CpG</H4>\n");
			pw.print("<P><IMG SRC=\"");
			pw.print(HtmlChartUrl(cycleCounters, Arrays.asList("CG"),0.0,1.0));
			pw.print("\">\n");
			pw.print("</P>\n");
		
		pw.print("<H4>CpH</H4>\n");
		pw.print("<P><IMG SRC=\"");
		pw.print(HtmlChartUrl(cycleCounters, Arrays.asList("CA","CY"),0.0,1.0));
		pw.print("\">\n");
		pw.print("</P>\n");
		
		pw.print("<H4>CpH zoom</H4>\n");
		pw.print("<P><IMG SRC=\"");
		//this.numCycles = 20;
		pw.print(HtmlChartUrl(cycleCounters, Arrays.asList("CA","CY"),0.96,1.0));
		pw.print("\">\n");
		pw.print("</P>\n");
		}
		else
		{
			pw.print("</P>\n");
			pw.print("<H4>% methylCs by cycle</H4>\n");
			pw.print("<P><IMG SRC=\"");
			pw.print(HtmlChartUrl(cycleCounters, null,0.005*(70/this.numCycles), 0.02*(70/this.numCycles)));
			pw.print("\">\n");
			pw.print("</P>\n");

		}

		pw.close();

//		// CSV output
//		System.out.println("\nCpG cycle counter");
//		System.out.println(cpgCycleCounter.toString());
//
//		System.out.println("\nCpH cycle counter");
//		System.out.println(cphCycleCounter.toString());

	}

	protected String HtmlChartUrl(Map<String,ByCycleCounter> counters, List<String> contextsToDisplay,
			double minVal, double maxVal)
	throws Exception
	{
		String out = null;

		List<Plot> plots = new ArrayList<Plot>(20);

		if (contextsToDisplay == null) contextsToDisplay = new ArrayList<String>(counters.keySet());
		
		try
		{
			// Get plots
			int contextCount = 1;
			int numCycles = 0;
			for (String context : contextsToDisplay)
			{
				ByCycleCounter counter = counters.get(context);
				numCycles = Math.max(numCycles, counter.numCycles);
				numCycles = Math.min(numCycles, this.numCycles);
				System.err.printf("numCycles(%d) = max(this.numcycles=%d, counter.numCycles=%d\n",numCycles, this.numCycles, counter.numCycles);

				// Once with filter, once without
				Plot noFilter = counter.toPlot(false, minVal,maxVal, this.relativeFreqs, numCycles);
				noFilter.setColor(GetContextColor(context, false));
				noFilter.setLegend(String.format("%s no filter", context));
				//noFilter.addShapeMarkers(Shape.CIRCLE, Color.BLACK, 3);
				plots.add(noFilter);

				if (!noFilters && !this.relativeFreqs)
				{	
					Plot withFilter = counter.toPlot(true, minVal,maxVal, this.relativeFreqs, numCycles);
					withFilter.setColor(GetContextColor(context, true));
					withFilter.setLegend(String.format("%s with filter", context));
					//withFilter.addShapeMarkers(Shape.SQUARE, Color.BLACK, 3);
					plots.add(withFilter);
				}

				contextCount++;

			}

			LineChart chart = GCharts.newLineChart(plots);
			chart.setSize(600, 200);
			chart.setDataEncoding(DataEncoding.EXTENDED);

			AxisLabels xAxis = AxisLabelsFactory.newNumericRangeAxisLabels(1, numCycles);
			chart.addXAxisLabels(xAxis);
			chart.addXAxisLabels(AxisLabelsFactory.newAxisLabels("Cycle number", 50.0));

			// This function seems to be buggy.
			//yAxis = AxisLabelsFactory.newNumericRangeAxisLabels(minVal, maxVal);

			AxisLabels yAxis;
			List<String> ylabels = new ArrayList<String>(10);
			double step = (maxVal-minVal)/5.0;
			for (double i = minVal; i <= maxVal; i = (Math.floor((i+step)*10000.0)/10000.0))
			{
				ylabels.add(String.format("%.2f", i));
				System.err.println("Adding label " + i);
			}
			yAxis = AxisLabelsFactory.newAxisLabels(ylabels);
			chart.addYAxisLabels(yAxis);
			chart.addRightAxisLabels(yAxis);
			chart.addYAxisLabels(AxisLabelsFactory.newAxisLabels("% conv", 50.0));
			

			chart.setTitle(ListUtils.excelLine(contextsToDisplay));
			out = chart.toURLString();
		}
		catch (Exception e)
		{
			System.err.println("Could not make chart: " + e.toString());
			e.printStackTrace(System.err);
		}


		return out;	
	}
	
	private static Color GetContextColor(String context, boolean filtered) {
		
		Color out = Color.ORANGE;
		if (context.equalsIgnoreCase("CG"))
		{
			out = (filtered) ? Color.GRAY : Color.BLACK;
		}
		else if (context.equalsIgnoreCase("CA"))
		{
			out = (filtered) ? Color.LIGHTGREEN : Color.GREEN ;
		}
		else if (context.equalsIgnoreCase("CT"))
		{
			out = (filtered) ?  Color.PINK : Color.RED;
		}
		else if (context.equalsIgnoreCase("CY"))
		{
			out = (filtered) ?  Color.PINK : Color.RED;
		}
		else if (context.equalsIgnoreCase("CC"))
		{
			out = (filtered) ?  Color.BEIGE : Color.BROWN;
		}
		
		return out;
	}

	static boolean isCytosine(int pos, String refStr, String seqStr)
	{
		char refC = refStr.charAt(pos);
		char seqC = seqStr.charAt(pos);
		
		return ((refC == 'C') && ((seqC == 'C') || (seqC == 'T'))); 
	}
	
	static boolean isCpg(int pos, String refStr, String seqStr)
	{
		if (pos >= (refStr.length()-1)) return false; // At the last character
		
		char refCnext = refStr.charAt(pos+1);
		char seqCnext = seqStr.charAt(pos+1);
		
		return ( isCytosine(pos,refStr,seqStr) && (refCnext == 'G') && (seqCnext == 'G') );
	}
	
	static String cytosineContext(int pos, String refStr, String seqStr)
	{
		if (pos >= (refStr.length()-1)) return "end"; // At the last character
		
		char refCcur = refStr.charAt(pos);
		char refCnext = refStr.charAt(pos+1);
		char seqCnext = seqStr.charAt(pos+1);
				
		// This is tricky.  We want to use the following base in the sequence,
		// unless it's a T.  If it's a T , you really can't tell what it is
		// unless you look at reads on the opposite strand.  And you don't want
		// to call CpCs either, because they might be biased for mCs.  So to just
		// show the information we're confident in , just show Y vs. A vs. G
		char next = seqCnext;
		if ((next == 'T') || (next == 'C')) next = 'Y';
//		if (next == 'T')
//		{
//			if (refCnext == 'C')
//			{
//				// You really don't know at this point which it is, unless
//				// you look at reads on the opposite strand
//				next = 'Y';
//			}
//		}

		return (String.format("%c%c",refCcur, next)).toUpperCase();
	}
	
	static boolean isConverted(int pos, String refStr, String seqStr)
	{
		char refC = refStr.charAt(pos);
		char seqC = seqStr.charAt(pos);
		boolean converted = ((refC == 'C') && (seqC == 'T'));
		//if (pos==0) System.err.printf("Pos %d, refC=%c, seqC=%c, converted=%s\n", pos, refC, seqC, converted);
		
		return converted;
	}

	static int boolToInd(boolean bool)
	{
		return ( bool ? 0 : 1 );
	}

	static String headerLine()
	{
		return "\tpre\t\t\t\tpost" + "\n" + "\tcpg\t\tcph\t\tcpg\t\tcph" + "\n" + "convS\tc\tnc\tc\tnc\tc\tnc\tc\tnc"; 
	}
	
	
	public class ByCycleCounter
	{
		// ind 1 = Pre or post.
		// ind 2 = conv/non-converted
		// ind 3 = cycle
		int[][][] counter; 
		int numCycles = 0;
		
		
		/**
		 * @param numCycles
		 */
		public ByCycleCounter(int numCycles) {
			this.numCycles = numCycles;
			counter = new int[2][2][numCycles];
		}

		void increment(boolean pre, boolean converted, int cycle)
		{
			//if (converted && (cycle==0)) System.err.printf("Incrementing cycle %d, pre=%s, converted=%s\n",cycle,pre,converted);
			if (cycle < numCycles)
			{
				counter[boolToInd(pre)][boolToInd(converted)][cycle]++;
			}
		}

		public Plot toPlot(boolean withFilter, double minScale, double maxScale, boolean relativeFreqs)
		throws Exception
		{
			return this.toPlot(withFilter, minScale, maxScale, relativeFreqs, -1);
		}
		
		public Plot toPlot(boolean withFilter, double minScale, double maxScale, boolean relativeFreqs, int maxCycles)
		throws Exception
		{
			System.err.printf("To plot with %f-%f\n",minScale,maxScale);
			Data data = this.toData(withFilter, minScale, maxScale, relativeFreqs, maxCycles);
			Plot plot = Plots.newPlot(data);
			

//			Color c;
//			switch (type)
//			{
//			case FW:
//				c = Color.BLUE; break;
//			case REV:
//				c = Color.RED; break;
//			default:
//				c = Color.CORNFLOWERBLUE; break;
//			}
//
//			// Bug in charts4j makes it so this doesn't work with 
//			// charts with lots of data points.  There is a proposed 
//			// workaround which i added to the project.
//			if (data.getSize() < 250)
//			{
//				plot.addShapeMarkers(Shape.DIAMOND, Color.BLACK, 4);
//			}
//			
//			// Add a line halfway
//			int midpoint = (int)Math.round((double)data.getSize()/2.0);
//			plot.addShapeMarker(Shape.VERTICAL_LINE_FULL, Color.BLACK,1,midpoint);
//			
//			plot.setColor(c);
			
			return plot;		
		}
		
		public Data toData(boolean withFilter, double minScale, double maxScale, boolean relativeFreqs)
		throws Exception
		{
			return this.toData(withFilter, minScale, maxScale, relativeFreqs, -1); 
		}

		public Data toData(boolean withFilter, double minScale, double maxScale, boolean relativeFreqs, int maxCycles)
		throws Exception
		{
			double[][][] dblCount = MatUtils.intMatToDouble(counter);
			
			double[][] convNonconv = (withFilter) ? dblCount[boolToInd(false)] : MatUtils.sumMats(dblCount[0], dblCount[1]);
			double[] data = MatUtils.divVects(convNonconv[boolToInd(true)], MatUtils.vectSum(convNonconv[0], convNonconv[1]));
			
			// Why was this in here?
			//data = MatUtils.vectSum(convNonconv[0], convNonconv[1]);
			
			if (relativeFreqs)
			{
				double[] meCounts = convNonconv[boolToInd(true)];
				double totalMe = MatUtils.nanSum(meCounts);
				data = MatUtils.divVect(meCounts, totalMe);
			}
			
			if (maxCycles>0)
			{
				data = Arrays.copyOfRange(data, 0, maxCycles);
			}
			
			// Nans not allowed
			MatUtils.nansToVal(data, minScale-1);
			System.err.println("Making data: " + ListUtils.excelLine(data));

			
			Data out = DataUtil.scaleWithinRange(minScale, maxScale, data);
			
			return out;
		}

		public String toString()
		{
			// Do it with a stringBuf so it's faster
			clearStrBuf();
			
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 4; j++)
				{
//					sb.append(String.format(", args)
//							"byCycle line: " + i + ", " + j + ",");
					for (int k = 0; k < numCycles; k++)
					{
						sb.append(counter[i][j][k] + ",");
					}
					sb.append("\n");
				}
			}			
			
			return sb.toString();
		}
	
	}
	
	
	
//	public class CytosineCounter
//	{
//		
//		// ind 1 = Pre or post.
//		// ind 2 = CpG or CpH:   0 = A, 1 = C, 2 = G, 3 = T 
//		// ind 3 = converted/non-converted
//		int[][][] counter = new int[2][4][2];
//	
//		
//		void increment(boolean pre, char base, boolean converted)
//		{
//			//System.err.println("Pre = " + pre);
////			counter[boolToInd(pre)][boolToInd(cpg)][boolToInd(converted)]++;
//			counter[boolToInd(pre)][base][boolToInd(converted)]++;
//		}
//		
//		public String toString()
//		{
//			// Do it with a stringBuf so it's faster
//			clearStrBuf();
//			
//			for (int i = 0; i < 2; i++)
//			{
//				for (int j = 0; j < 2; j++)
//				{
//					for (int k = 0; k < 2; k++)
//					{
//						sb.append(counter[i][j][k] + "\t");
//					}
//				}
//			}
//			
//			return sb.toString();
//		}
//	
//	}
	
	public static int indToNuc(int ind)
	{
		char nuc = 0;
		switch (ind)
		{
		case 0:
			nuc = 'A'; break;
		case 1:
			nuc = 'C'; break;
		case 2:
			nuc = 'T'; break;
		case 3:
			nuc = 'G'; break;
		case 4:
			nuc = 'N'; break;
		default:
			System.err.printf("nucToInd: %d is not a valid nucleotide index", ind);
			System.exit(1);
		}
		return nuc;
	}
	
	public static int nucToInd(char nuc)
	{
		int out = -1;
		switch (nuc)
		{
		case 'a':
		case 'A':
			out = 0; break;
		case 'c':
		case 'C':
			out = 1; break;
		case 't':
		case 'T':
			out = 2; break;
		case 'G':
		case 'g':
			out = 3; break;
		case 'N':
		case 'n':
			out = 4; break;
		default:
			System.err.printf("nucToInd: %c is not a valid nucleotide", nuc);
			System.exit(1);
		}
		return out;
	}
	
	
}
