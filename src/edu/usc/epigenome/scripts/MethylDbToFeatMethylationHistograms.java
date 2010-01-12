package edu.usc.epigenome.scripts;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.seq.StrandedFeature;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.usckeck.genome.ChromFeatures;
import org.usckeck.genome.GFFUtils;

import com.googlecode.charts4j.Color;


import edu.usc.epigenome.genomeLibs.Charts4jUtils;
import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAligner;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerEachfeat;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRangeWithRefpoint;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgSummarizer;


public class MethylDbToFeatMethylationHistograms {

	private static List<Color> colors = Arrays.asList(Color.BLUE, Color.RED, Color.GREEN, Color.ORANGE, Color.BEIGE, Color.BLACK, Color.AQUAMARINE);

	private static final String C_USAGE = "Use: MethylDbToFeatMethylationHistograms -minCTreads 0 -maxOppStrandAfrac 0.1 -maxFeatSize 10 -flankSize 2000 " +
	"-outputPrefix outputTag sample1_tablePrefix sample2_tablePrefix ... , feats1.gtf feats2.gtf ...";
	
	@Option(name="-featFilter",usage="We will take the intersection with this feature. Must be a featType in the features table")
	protected List<String> featFilters = new ArrayList<String>(5);
	@Option(name="-maxFeatSize",usage="maximum size of features to include (default Inf)")
    protected int maxFeatSize = Integer.MAX_VALUE;
    @Option(name="-flankSize",usage="bp flanking each side of the feature center (default 2000)")
    protected int flankSize = 2000;
    @Option(name="-outputPrefix",usage="Prefix for output files (default methylDb)")
    protected String outputPrefix = "methylDb";
    @Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
    protected int minCTreads = 0;
    @Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
    		"This sets a maximum number of observed As on the opposite strand (default 0.1)")
    protected double maxOppStrandAfrac = 0.1;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();



	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception	
	{
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
		new MethylDbToFeatMethylationHistograms().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception
	{
		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);
		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			System.err.println(C_USAGE);
			// print the list of available options
			parser.printUsage(System.err);
			System.err.println();
			return;
		}

		List<String> tablePrefixes = new ArrayList<String>(3);
		List<String> featFns = new ArrayList<String>(100);
		boolean onFeats = false;
		for (String arg : arguments)
		{
			if (arg.equals(","))
			{
				//System.err.println("Found comma!");
				onFeats = true;
			}
			else
			{
				List<String> l = (onFeats) ? featFns : tablePrefixes;
				l.add(arg);
			}
		}
		
		if( (tablePrefixes.size()<1) || (featFns.size()<1))
		{
			System.err.println(C_USAGE);
			parser.printUsage(System.err);
			return;
		}
		
		int nS = tablePrefixes.size();
		int nFeatTypes = featFns.size();

		// Setup writer
		PrintWriter writer = new PrintWriter(new FileOutputStream(String.format("%s.histograms.html", outputPrefix)));
		writer.println("<P>" + ListUtils.excelLine(args) + "</P>");

		// Go through summarizer types
		// Output matrices
		List<String> summarizerTypes = Arrays.asList("edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer",
				"edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgCoverageSummarizer");
//				"edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgDeaminationSummarizer");
//				"edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgNonconversionSummarizer");
	
		// Start the actual work.  Go through each GFF file.
		int onFeatType = 0;
		for (String featFn : featFns)
		{
			try
			{
			onFeatType++;
			ChromFeatures feats = new ChromFeatures(featFn, true);
			System.err.println("About to filter regions by size");
			feats = feats.filterBySize(0, this.maxFeatSize);


			String featsFnBase = (new File(featFn)).getName();
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info("On feature " + onFeatType + "/" + nFeatTypes + " " + featsFnBase);
			// The main work
			double[][][] methLevels = getMethLevels(summarizerTypes, feats, tablePrefixes);
			int nFeats = methLevels[0][0].length;
			writer.printf("<H2>%s (%d feats)</H2>\n<P></P>\n", featsFnBase, nFeats);


			for (int summ=0; summ < summarizerTypes.size() ; summ++  )
			{
				double min = 0.0;
				double max = 1.0;
				String summName = summarizerTypes.get(summ);
				if (summName.contains("Coverage"))
				{
					min = 0.0;
					max = 20.0;
				}
				else if (summName.contains("eamination") || summName.contains("onversion"))
				{
					min = 0.0;
					max = 1.0;
				}

				writer.append("<H4>");
				writer.append(summName);
				writer.append("</H4>\n");
				writer.append("<P>");
				for (int s = 0; s < nS; s++)
				{

					String tableName = tablePrefixes.get(s);
					tableName = tableName.replaceAll("methylCGsRich_", "");
					tableName = tableName.replaceAll("_", "");
					double[] means = methLevels[s][summ];
					String chartUrl = Charts4jUtils.histogramBarPlot(means, min, max, 10, 
							String.format("%s",tableName), colors.get(summ));
					writer.append(chartUrl);
				}
				writer.append("</P>\n");
			}
			}
			catch (Exception e)
			{
				System.err.println("Couldn't complete feature " + featFn + ":\n" + e.toString());
				e.printStackTrace();
			}

		}
		writer.close();
	}

	
	/**
	 * @param chrStr
	 * @param feats
	 * @param tablePrefixes
	 * @return double[i][j][k] where i is the table (sample), k is the summarizer type and j is the feature mean
	 * @throws Exception
	 */
	protected double[][][] getMethLevels(List<String> summTypes, ChromFeatures feats, List<String> tablePrefixes)
	throws Exception
	{
		// Setup the DB queries
		int nS = tablePrefixes.size();
		int nF = feats.num_features();
		int nSumm = summTypes.size();
		double[][][] out = new double[nS][nSumm][nF];
		MatUtils.initMat(out, Double.NaN);
		
		int[] featNum = new int[nSumm];

		for (String chrStr : MethylDbUtils.CHROMS)
		{
			int chr = feats.chrom_from_public_str(chrStr); 
			Iterator featit = feats.featureIterator(chr);
			System.err.println("Processing " + chrStr);
			FEAT: while (featit.hasNext())
			{

				SimpleGFFRecord rec = (SimpleGFFRecord)featit.next();
				StrandedFeature.Strand featStrand = rec.getStrand();
				String featName = null; // rec.getSeqName();
				int featS = rec.getStart();
				int featE = rec.getEnd();
				int flankStart = featS - this.flankSize;
				int flankEnd = featE + this.flankSize;


				Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(String.format(
						"Fetching coords: chr=%s\tfeatS=%d\tfeatE=%d\tfeatStrand=%s\tflankS=%d\tflankEnd=%d\t\n",
						chrStr, featS, featE, ""+featStrand, flankStart, flankEnd));


				try 
				{
					// Setup query params
					MethylDbQuerier params = new MethylDbQuerier();
					for (String featFilter : this.featFilters)
					{
						params.addFeatFilter(featFilter,flankSize);
					}
					params.setMinCTreads(this.minCTreads);
					params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
					params.addRangeFilter(chrStr,flankStart,flankEnd);

					// Setup summarizers
					for (int summ = 0 ; summ < nSumm; summ++)
					{
						Class<CpgSummarizer> summClass = (Class<CpgSummarizer>) Class.forName(summTypes.get(summ));
						CpgSummarizer[] summarizers = new CpgSummarizer[nS];
						for (int i=0; i<nS;i++)
						{
							summarizers[i] = summClass.newInstance();
							summarizers[i].init(params);
						}
						
						// And stream through summarizers
						CpgIteratorMultisample cpgit = new CpgIteratorMultisample(params, tablePrefixes);
						while (cpgit.hasNext()) 
						{
							Cpg[] cpgs = cpgit.next();

							for (int i = 0; i < nS; i++)
							{
								summarizers[i].streamCpg(cpgs[i]);
							}
						}

						// And set output
						for (int i = 0; i < nS; i++) out[i][summ][featNum[summ]] = summarizers[i].getValMean();
						featNum[summ]++;
					}

				}
				catch (Exception e)
				{
					Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe(String.format(
							"PROBLEM WITH FEATURE %s coords: chr=%s\tfeatS=%d\tfeatE=%d\tfeatStrand=%s\tflankS=%d\tflankEnd=%d\n%s\n",
							featName, chrStr, featS, featE, ""+featStrand, flankStart, flankEnd, e.toString()));
					e.printStackTrace();
					System.exit(1);

				}



			}
		}

		// We might not have done all the chromosomes, so downsize the array
		for (int i = 0; i < nS; i++)
		{
			for (int j = 0; j < nSumm; j++)
			{
				out[i][j] = Arrays.copyOf(out[i][j], featNum[0]);
			}
		}
		
		return out;

	}

	
	
}
