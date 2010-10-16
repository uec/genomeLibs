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

import com.mallardsoft.tuple.Pair;
import com.mallardsoft.tuple.Tuple;


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
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalker;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerParams;


public class MethylDbToWindowMeth {

	private static final String C_USAGE = "Use: MethylDbToWindowMeth sample1_tablePrefix sample2_tablePrefix ... , feats1.gtf feats2.gtf ...";
	
	@Option(name="-maxFeatSize",usage="maximum size of features to include (default Inf)")
    protected int maxFeatSize = Integer.MAX_VALUE;
	@Option(name="-useSpatialCpgWeighting",usage="If set, weights a CpG by the spatial distance it covers, otherwise weight all CpGs equally (default false)")
    protected boolean useSpatialCpgWeighting = false;
	@Option(name="-windSize",usage="add one or more window sizes to be tested in addition to the domain of the feature itself")
    protected List<Integer> windSizes = new ArrayList<Integer>();
    @Option(name="-outputPrefix",usage="Prefix for output files (default methylDb)")
    protected String outputPrefix = "methylDb";
    @Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
    protected int minCTreads = 0;
    @Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
    		"This sets a maximum number of observed As on the opposite strand (default 0.1)")
    protected double maxOppStrandAfrac = 0.1;
	@Option(name="-featFilter",usage="We will take the intersection with this feature (only for -windSize windows, NOT for primary feats). Must be a featType in the features table")
	protected List<String> featFilters = new ArrayList<String>(5);
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();


	// class vars
	int fCurFeatInd = 0;

	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception	
	{
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
		new MethylDbToWindowMeth().doMain(args);
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
		
		windSizes.add(0, new Integer(0));

		
		int nS = tablePrefixes.size();
		int nFeatTypes = featFns.size();

		// Setup writer
		PrintWriter writer = new PrintWriter(new FileOutputStream(String.format("%s.windAvs.csv", outputPrefix)));
		

	
		// Start the actual work.  Go through each GFF file.
		int onFeatType = 0;
		for (String featFn : featFns)
		{
			onFeatType++;
			ChromFeatures feats = new ChromFeatures(featFn, true);
			System.err.println("About to filter regions by size");
			feats = feats.filterBySize(0, this.maxFeatSize);
			int nFeats = feats.num_features(); 
			
//			PrintWriter sortWriter = new PrintWriter(new FileOutputStream(
//					String.format("%s.flank%d.featType%d.sortVals.csv", outputPrefix, this.flankSize, onFeatType)));
//			

			for (String chrStr : MethylDbUtils.CHROMS) // Arrays.asList("chr11")) //
			{
				processChrom(chrStr, feats, tablePrefixes, writer);
			}


		}
		writer.close();
	}


	protected void processChrom(String chrStr, ChromFeatures feats, List<String> tablePrefixes, PrintWriter writer)
	throws Exception
	{
		// Setup the DB queries
		int chr = feats.chrom_from_public_str(chrStr); 
		int nS = tablePrefixes.size();
		
		// Make a CpgMethLevelSummarizer for each table
		CpgMethLevelSummarizer[] summarizers = new CpgMethLevelSummarizer[nS];
		for (int i =0; i<nS; i++)
		{
			summarizers[i] = new CpgMethLevelSummarizer();
		}
	

		Iterator featit = feats.featureIterator(chr);
		System.err.println("Processing " + chrStr);
		int featNum = 0;
		FEAT: while (featit.hasNext())
		{
			featNum++;

			SimpleGFFRecord rec = (SimpleGFFRecord)featit.next();
			StrandedFeature.Strand featStrand = rec.getStrand();
			String featName = null; // rec.getSeqName();
			int featS = rec.getStart();
			int featE = rec.getEnd();

			for (Integer windSize : windSizes)
			{

				try 
				{

					int featS_w = featS;
					int featE_w = featE;
					if (windSize>0)
					{
						double windFlank = (double)windSize / 2.0;
						double windC = Math.round(MatUtils.mean(featS, featE));
						featS_w = (int)Math.round(windC - windFlank);
						featE_w = (int)Math.round(windC + windFlank);
					}
					
					


					// Meth
					MethylDbQuerier params = new MethylDbQuerier();

					// Only use features for window size > 0
					if (windSize>0)
					{
						for (String featFilter : this.featFilters)
						{
							params.addFeatFilter(featFilter,0);
						}
					}
					params.setMinCTreads(this.minCTreads);
					params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
					params.addRangeFilter(chrStr,featS_w,featE_w);

					CpgIteratorMultisample cpgit = new CpgIteratorMultisample(params, tablePrefixes);
					
					// Reset summarizers
					for (int i =0; i<nS; i++)
					{
						summarizers[i].reset();
					}


					int nCpgs = 0;
					while (cpgit.hasNext()) 
					{
						Cpg[] cpgs = cpgit.next();

						int chromPos = cpgs[0].chromPos;
						StrandedFeature.Strand cpgStrand = cpgs[0].getStrand();


						for (int i=0; i<nS; i++)
						{
							summarizers[i].streamCpg(cpgs[i]);
						}
						
						//						System.err.println(String.format("\tON CPG %s", cpgs[0].toString()));
						nCpgs++;
					}


					if ((featNum % 100)==0)
					{
						System.err.print(String.format(
								"ON WINDSIZE=%d, FEATURE %d coords: chr=%s\tfeatS=%d\tfeatE=%d\tflankS=%d\tflankE=%d\tfeatStrand=%s\n\t%d Cpgs\n",
								windSize, featNum, chrStr, featS, featE, featS_w, featE_w, ""+featStrand, nCpgs));
					}
					
					for (int i =0; i<nS; i++)
					{
						writer.printf("%.2f,", summarizers[i].getValMean(this.useSpatialCpgWeighting));
					}


				}
				catch (Exception e)
				{
					Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe(String.format(
							"PROBLEM WITH FEATURE %d coords: chr=%s\tfeatS=%d\tfeatE=%d\tfeatStrand=%s\n%s\n",
							featNum, chrStr, featS, featE, ""+featStrand,e.toString()));
					e.printStackTrace();
					System.exit(1);

				}

			}
			// Increment feat ind
			this.fCurFeatInd++;
			
			// New line
			writer.println();

		}

		

	}
	
	
	
}
