package edu.usc.epigenome.scripts;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.seq.StrandedFeature;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.usckeck.genome.ChromFeatures;


import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerEachfeat;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


public class MethylDbToMultisampleFeatAlignmentsStratified {

	private static final String C_USAGE = "Use: MethylDbToMultisampleFeatAlignmentsStratified -censor -alignToStart  -maxFeatSize 10 -skipUnoriented -flankSize 2000 " +
	"-outputPrefix outputTag -downsampleCols 500 -heatmap sample1_tablePrefix sample2_tablePrefix ... , feats1.gtf feats2.gtf ...";
	
	@Option(name="-skipUnoriented",usage="If set, skip any unoriented feature (default false)")
	protected boolean skipUnoriented = false;
	@Option(name="-combineStrands",usage="If set, combine strands into a single line")
	protected boolean combineStrands = false;
	@Option(name="-heatmap",usage="If set, make a java heatmap")
	protected boolean heatmap = false;
//	@Option(name="-noDeltas",usage="If set, do not output any delta plots")
//	protected boolean noDeltas = false;
	@Option(name="-censor",usage="If set, do not include points within the flank region but inside the feature region")
	protected boolean censor = false;
	@Option(name="-alignToStart",usage="If set, align to the left end (or 5' if available) of the feature.  Default is to align to center")
	protected boolean alignToStart = false;
	@Option(name="-alignToEnd",usage="If set, align to the right end (or 3' if available) of the feature.  Default is to align to center")
	protected boolean alignToEnd = false;
	@Option(name="-maxFeatSize",usage="maximum size of features to include (default Inf)")
    protected int maxFeatSize = Integer.MAX_VALUE;
    @Option(name="-flankSize",usage="bp flanking each side of the feature center (default 2000)")
    protected int flankSize = 2000;
    @Option(name="-downscaleCols",usage="Number of cols in ouput (downscaled from flank width)")
    protected int downscaleCols = 500;
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


	// class vars
	FeatAlignerEachfeat[][] fStatMats = null; // n-> sample number, m -> 0=readCount, 1=nCpGs, 2=mLevel
//	FeatAligner[] fDeltaMats = null; // n -> (nS*(nS-1))/2  (all pairwise)
//	FeatAligner fVarianceMat = null;
	int fCurFeatInd = 0;


	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception	
	{
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.INFO);
		new MethylDbToMultisampleFeatAlignmentsStratified().doMain(args);
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
		if (this.censor && !this.skipUnoriented)
		{
			System.err.println("With censoring, you should add -skipUnoriented");
			parser.printUsage(System.err);
			return;
		}
		
		int nS = tablePrefixes.size();
		int nFeatTypes = featFns.size();

		// Setup writer
		PrintWriter writer = new PrintWriter(new FileOutputStream(String.format("%s.charts.html", outputPrefix)));
		writer.println("<P>" + parser.toString() + "</P>");

	
		// Start the actual work.  Go through each GFF file.
		int onFeatType = 0;
		for (String featFn : featFns)
		{
			onFeatType++;
			ChromFeatures feats = new ChromFeatures(featFn, true);
			System.err.println("About to filter regions by size");
			feats = feats.filterBySize(0, this.maxFeatSize);
			int nFeats = feats.num_features(); 
			
			
			
			// Create arrays
			System.err.println("About to initialize aligners");
			fStatMats = new FeatAlignerEachfeat[nS][3];
//			fDeltaMats = new FeatAligner[(nS*(nS-1))/2];
//			fVarianceMat = new FeatAlignerAveraging(flankSize,false);
			int onDeltaMat = 0;
			for (int i = 0; i < nS; i++)
			{
				// 0=readCount, 1=nCpGs, 2=mLevel
//				fStatMats[i][0] = new FeatAlignerEachfeat(flankSize,false, nFeats,500);
//				fStatMats[i][1] = new FeatAlignerEachfeat(flankSize,true, nFeats,500);
				fStatMats[i][2] = new FeatAlignerEachfeat(flankSize,false, nFeats,this.downscaleCols);
//				for (int j = (i+1); j < nS; j++)
//				{
//					fDeltaMats[onDeltaMat++] = new FeatAlignerAveraging(flankSize,false);
//				}
			}
	
			for (String chrStr : MethylDbUtils.CHROMS)
			{
				processChrom(chrStr, feats, tablePrefixes, skipUnoriented);
			}

		
			// Output matrices
			String featsFnBase = (new File(featFn)).getName();
			writer.printf("<H1>%s</H1>\n<P></P>\n", featsFnBase);
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info("On feature " + onFeatType + "/" + nFeatTypes + " " + featsFnBase);
		

			// M-levels
			for (int i = 0; i < nS; i++)
			{
				String tablePrefix = tablePrefixes.get(i);
				writer.printf("<H4>%s</H4>\n", tablePrefix);
				writer.println(this.fStatMats[i][2].htmlChart(!this.combineStrands, true, true));
				
				double[] colorMinMax = {0.0,1.0};
				if (this.heatmap) this.fStatMats[i][2].launchSwingHeatmap(colorMinMax);
			}

//			if (!this.noDeltas)
//			{
//				// Pairwise deltas
//				int onMat = 0;
//				for (int i = 0; i < nS; i++)
//				{
//					String aPrefix = tablePrefixes.get(i);
//					for (int j = (i+1); j < nS; j++)
//					{
//						String bPrefix = tablePrefixes.get(j);
//
//						writer.printf("<H4>Abs diff %s vs. %s</H4>\n", aPrefix, bPrefix);
//						writer.println(this.fDeltaMats[onMat++].htmlChart(!this.combineStrands, true, false));
//					}
//				}
//
//				// Variance
//				writer.printf("<H4>Variance</H4>\n");
//				writer.println(this.fVarianceMat.htmlChart(!this.combineStrands, true, false));
//			}

			


		}
		writer.close();
	}
	
	
	protected void processChrom(String chrStr, ChromFeatures feats, List<String> tablePrefixes, boolean skipUnoriented)
	throws Exception
	{
		// Setup the DB queries
		int chr = feats.chrom_from_public_str(chrStr); 
		int nS = tablePrefixes.size();

		Iterator featit = feats.featureIterator(chr);
		System.err.println("Processing " + chrStr);
		int featNum = 0;
		while (featit.hasNext())
		{
			featNum++;
			
			GFFRecord rec = (GFFRecord)featit.next();
			StrandedFeature.Strand featStrand = rec.getStrand();
			String featName = null; // rec.getSeqName();

			if (skipUnoriented)
			{
				// Don't use those without orientation
				if (featStrand == StrandedFeature.UNKNOWN) continue;
			}
			else
			{
				if (featStrand == StrandedFeature.UNKNOWN) featStrand = StrandedFeature.POSITIVE;
			}
				
			// Do we want to center align?  Do we want to censor
			int featS = rec.getStart();
			int featE = rec.getEnd();
			int featCenter = (int)Math.round(((double)featE+(double)featS) / 2.0);

			// Get the alignment point.
			int alignmentPoint;
			if (this.alignToStart)
			{
				alignmentPoint = (featStrand == StrandedFeature.POSITIVE) ? featS : featE;
			}
			else if (this.alignToEnd)
			{
				alignmentPoint = (featStrand == StrandedFeature.POSITIVE) ? featE : featS;
			}
			else
			{
				alignmentPoint = featCenter;
			}
			
			// And the flank endpoints depends on censoring
			int flankStart, flankEnd;
			if (!this.censor)
			{
				flankStart = alignmentPoint-flankSize;
				flankEnd = alignmentPoint+flankSize;
			}
			else
			{
				flankStart = Math.max(featS, alignmentPoint-flankSize);
				flankEnd = Math.min(featE, alignmentPoint+flankSize);
				
				// Censoring is relative to feature strand ONLY if we aligne to 
				// start.  If we align to center, we censor on both sides.
				if (this.alignToStart)
				{
					if (featStrand == StrandedFeature.NEGATIVE)
					{
						flankEnd = alignmentPoint + flankSize;  // Don't censor 5'
					}
					else
					{
						flankStart = alignmentPoint - flankSize; // Don't censor 3'
					}
				}
				else if (this.alignToEnd)
				{
					if (featStrand == StrandedFeature.NEGATIVE)
					{
						flankStart = alignmentPoint - flankSize; // Don't censor 3'
					}
					else
					{
						flankEnd = alignmentPoint + flankSize;  // Don't censor 5'
					}
				}
			}
			
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(String.format(
					"Fetching coords: censor=%s\talignToStart=%s\tchr=%s\tfeatS=%d\tfeatE=%d\tfeatStrand=%s\talignmentPoint=%d\tflankS=%d\tflankEnd=%d\t\n",
					""+this.censor, ""+this.alignToStart, chrStr, featS, featE, ""+featStrand, alignmentPoint, flankStart, flankEnd));
					
	
			try 
			{

			// Meth
			MethylDbQuerier params = new MethylDbQuerier();
			params.setMinCTreads(this.minCTreads);
			params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
			params.addRangeFilter(chrStr,flankStart,flankEnd);
			CpgIteratorMultisample cpgit = new CpgIteratorMultisample(params, tablePrefixes);
			while (cpgit.hasNext()) 
			{
				Cpg[] cpgs = cpgit.next();
				
				int chromPos = cpgs[0].chromPos;
				StrandedFeature.Strand cpgStrand = cpgs[0].getStrand();
				
				
				// First individual stats
				double mLevels[] = new double[nS];
				for (int i = 0; i < nS; i++)
				{
					double mLevel = cpgs[i].fracMeth(params.getUseNonconversionFilter());
					mLevels[i] = mLevel;
					this.fStatMats[i][2].addAlignmentPos(
							chromPos,
							(cpgStrand == StrandedFeature.NEGATIVE) ? Double.NaN : mLevel,
									(cpgStrand == StrandedFeature.NEGATIVE) ? mLevel: Double.NaN,
											featName, chrStr, alignmentPoint, featStrand);
				}

//				if (!this.noDeltas)
//				{
//					// Then variance
//					double mVar = MatUtils.nanVariance(mLevels);
//					this.fVarianceMat.addAlignmentPos(
//							chromPos,
//							(cpgStrand == StrandedFeature.NEGATIVE) ? Double.NaN : mVar,
//									(cpgStrand == StrandedFeature.NEGATIVE) ? mVar: Double.NaN,
//											featName, chrStr, alignmentPoint, featStrand);
//
//					// Then pairwise.
//					int onMat = 0;
//					for (int i = 0; i < nS; i++)
//					{
//						for (int j = (i+1); j < nS; j++)
//						{
//							double absDiff = Math.abs(mLevels[i]-mLevels[j]);
//							this.fDeltaMats[onMat++].addAlignmentPos(
//									chromPos,
//									(cpgStrand == StrandedFeature.NEGATIVE) ? Double.NaN : absDiff,
//											(cpgStrand == StrandedFeature.NEGATIVE) ? absDiff : Double.NaN,
//													featName, chrStr, alignmentPoint, featStrand);
//						}
//					}
//				}

			}
			}
			catch (Exception e)
			{
				Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe(String.format(
						"PROBLEM WITH FEATURE %d coords: censor=%s\talignToStart=%s\tchr=%s\tfeatS=%d\tfeatE=%d\tfeatStrand=%s\talignmentPoint=%d\tflankS=%d\tflankEnd=%d\n%s\n",
						featNum, ""+this.censor, ""+this.alignToStart, chrStr, featS, featE, ""+featStrand, alignmentPoint, flankStart, flankEnd, e.toString()));
				e.printStackTrace();
				System.exit(1);
				
			}

			
			// Increment feat ind
			this.fCurFeatInd++;
		}
		
		

	}
	
	
	
}