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
import edu.usc.epigenome.genomeLibs.FeatDb.FeatDbQuerier;
import edu.usc.epigenome.genomeLibs.FeatDb.FeatIterator;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRangeWithRefpoint;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIteratorMultisample;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbUtils;


public class MethylDbToMultisampleFeatAlignmentsStratified {

	private static final String C_USAGE = "Use: MethylDbToMultisampleFeatAlignmentsStratified -censor -alignToStart  -maxFeatSize 10 -skipUnoriented -flankSize 2000 " +
	"-outputPrefix outputTag -downsampleCols 500 -heatmap -sortByExpression 'T14838A-N14838B' sample1_tablePrefix sample2_tablePrefix ... , feats1.gtf feats2.gtf ...";
	
	@Option(name="-skipUnoriented",usage="If set, skip any unoriented feature (default false)")
	protected boolean skipUnoriented = false;
	@Option(name="-combineStrands",usage="If set, combine strands into a single line")
	protected boolean combineStrands = false;
	@Option(name="-extendRead",usage="If set, extends the start position of the read by the indicated amount (in a strand-aware way)")
	protected int extendRead = 0;
	@Option(name="-heatmap",usage="If set, make a java heatmap")
	protected boolean heatmap = false;
	@Option(name="-nomatlab",usage="If set, don't make matlab output files (default false)")
	protected boolean nomatlab = false;
	@Option(name="-nometh",usage="If set, don't make meth files (default false)")
	protected boolean nometh = false;
	@Option(name="-nohtml",usage="If set, don't make html files (default false)")
	protected boolean nohtml = false;
	@Option(name="-nosort",usage="If set, don't sort at all (nullifies stratification). Default false")
	protected boolean nosort = false;
//	@Option(name="-noDeltas",usage="If set, do not output any delta plots")
//	protected boolean noDeltas = false;
	@Option(name="-featFilter",usage="We will take the intersection with this feature. Must be a featType in the features table")
	protected List<String> featFilters = new ArrayList<String>(5);
	@Option(name="-censor",usage="If set, do not include points within the flank region but inside the feature region")
	protected boolean censor = false;
	@Option(name="-alignToStart",usage="If set, align to the left end (or 5' if available) of the feature.  Default is to align to center")
	protected boolean alignToStart = false;
	@Option(name="-featTableCenters",usage="Use counts from the feature table rather than methyl table. In this case, sample1_tab prefix is actually a featType")
	protected boolean featTableCenters = false;
	@Option(name="-featTableStarts",usage="Use counts from the feature table rather than methyl table. In this case, sample1_tab prefix is actually a featType")
	protected boolean featTableStarts = false;
	@Option(name="-featTableStartsAndEnds",usage="Use counts from the feature table rather than methyl table. In this case, sample1_tab prefix is actually a featType")
	protected boolean featTableStartsAndEnds = false;
	@Option(name="-featTableCompleteExtents",usage="Use counts from the feature table rather than methyl table. In this case, sample1_tab prefix is actually a featType")
	protected boolean featTableCompleteExtents = false;
	@Option(name="-readCounts",usage="Output read counts")
	protected boolean readCounts = false;
	@Option(name="-readCountsIgnoreTotalReads",usage="Used only if -readCounts is set. Instead of using totalReads field, counts each database row as 1 (useful for CpG densities). Default false")
	protected boolean readCountsIgnoreTotalReads = false;
	@Option(name="-alignToEnd",usage="If set, align to the right end (or 3' if available) of the feature.  Default is to align to center")
	protected boolean alignToEnd = false;
	@Option(name="-includeCoords",usage="If set, include coords in main csv file (default false)")
	protected boolean includeCoords = false;
	@Option(name="-maxFeatSize",usage="maximum size of features to include (default Inf)")
    protected int maxFeatSize = Integer.MAX_VALUE;
    @Option(name="-flankSize",usage="bp flanking each side of the feature center (default 2000)")
    protected int flankSize = 2000;
    @Option(name="-downscaleCols",usage="Number of cols in ouput (downscaled from flank width)")
    protected int downscaleCols = 2000;
    @Option(name="-outputPrefix",usage="Prefix for output files (default methylDb)")
    protected String outputPrefix = "methylDb";
    @Option(name="-sortByExpression",usage="An expression from the gene expression table, used to sort all features that can be linked by RefSeq")
    protected String sortByExpression = null;
    @Option(name="-minCTreads",usage="Minimum number of C or T reads to count as a methylation value")
    protected int minCTreads = 0;
    @Option(name="-maxOppStrandAfrac",usage="As on the opposite strand are evidence for mutation or SNP. " +
    		"This sets a maximum number of observed As on the opposite strand (default 0.1)")
    protected double maxOppStrandAfrac = 0.1;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();


	// class vars
	FeatAlignerEachfeat[][] fStatMats = null; // n-> sample number, m -> 0=readCount, 1=nCpGs, 2=mLevel, 3=cpg weights
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
		
		if (this.includeCoords && !this.nosort)
		{
			System.err.println("If -includeCoords is specified, also must specify -nosort");
			System.err.println(C_USAGE);
			parser.printUsage(System.err);
			return;
			
		}
		
		if (this.alignToEnd && this.alignToEnd)
		{
			System.err.println("Can not include both -alignToStart and -alignToEnd");
			System.err.println(C_USAGE);
			parser.printUsage(System.err);
			return;
			
		}

		if ( ((this.featTableCenters?1:0) +
				(this.featTableCompleteExtents?1:0) + 
				(this.featTableStarts?1:0) + 
				(this.featTableStartsAndEnds?1:0)) > 1) 
		{
			System.err.println("Only one of -featTableStarts, -featTableStartsAndEnds, -featTableCenters and -featTableCompleteExtents can be specified");
			System.err.println(C_USAGE);
			parser.printUsage(System.err);
			return;
		}
		else if (this.featTableCenters || this.featTableCompleteExtents || this.featTableStarts || this.featTableStartsAndEnds)
		{
			this.readCounts = true;
			this.nometh = true;
		}
		

		// I don't think this is necessary.
		//		if (this.censor && !this.skipUnoriented)
//		{
//			System.err.println("With censoring, you should add -skipUnoriented");
//			parser.printUsage(System.err);
//			return;
//		}
		
		int nS = tablePrefixes.size();
		int nFeatTypes = featFns.size();

		// Setup writer
		PrintWriter writer = new PrintWriter(new FileOutputStream(String.format("%s.flank%d.charts.html", outputPrefix, this.flankSize)));
		writer.println("<P>" + ListUtils.excelLine(args) + "</P>");

	
		// Start the actual work.  Go through each GFF file.
		int onFeatType = 0;
		for (String featFn : featFns)
		{
			onFeatType++;
			ChromFeatures feats = new ChromFeatures(featFn, true);
			System.err.printf("About to filter regions by size, starting with %d\n", feats.num_features());
			feats = feats.filterBySize(0, this.maxFeatSize);
			System.err.printf("Filter regions by size, found %d between %d and %d\n", feats.num_features(), 0, this.maxFeatSize);
			int nFeats = feats.num_features(); 
			
			PrintWriter sortWriter = new PrintWriter(new FileOutputStream(
					String.format("%s.flank%d.featType%d.sortVals.csv", outputPrefix, this.flankSize, onFeatType)));
			
			// Create arrays
			System.err.println("About to initialize aligners");
			fStatMats = new FeatAlignerEachfeat[nS][4];
//			fDeltaMats = new FeatAligner[(nS*(nS-1))/2];
//			fVarianceMat = new FeatAlignerAveraging(flankSize,false);
			int onDeltaMat = 0;
			for (int i = 0; i < nS; i++)
			{
				// 0=readCount, 1=nCpGs, 2=mLevel, 3= CpG weights
				if (readCounts) fStatMats[i][0] = new FeatAlignerEachfeat(flankSize,!this.censor, nFeats,this.downscaleCols); // Changed this to start with NaN.  Now we explicitly zero out features to allow censoring
//				fStatMats[i][1] = new FeatAlignerEachfeat(flankSize,true, nFeats,500);
				if (!nometh) fStatMats[i][2] = new FeatAlignerEachfeat(flankSize,false, nFeats,this.downscaleCols);
				if (!nometh) fStatMats[i][3] = new FeatAlignerEachfeat(flankSize,false, nFeats,this.downscaleCols);
//				for (int j = (i+1); j < nS; j++)
//				{
//					fDeltaMats[onDeltaMat++] = new FeatAlignerAveraging(flankSize,false);
//				}
			}
	
			for (String chrStr : MethylDbUtils.CHROMS) // Arrays.asList("chr11")) //
			{
				processChrom(chrStr, feats, tablePrefixes, skipUnoriented);
			}

		
			// Output matrices
			String featsFnBase = (new File(featFn)).getName();
			writer.printf("<H1>%s</H1>\n<P></P>\n", featsFnBase);
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info("On feature " + onFeatType + "/" + nFeatTypes + " " + featsFnBase);
		

			// M-levels
			Double[] sortVals = null;
			int ind = (this.readCounts) ? 0 : 2;
			for (int i = 0; i < nS; i++)
			{
				if (!this.nosort)
				{
					if (this.sortByExpression != null)
					{
							sortVals = this.fStatMats[i][ind].sortRowsBySortVals();
					}
					else if (sortVals == null)
					{
						int nCols = this.fStatMats[i][ind].numCols();
						int colsStart = 0;
						int colsEnd = nCols-1;
						if (this.censor)
						{
							int midPoint = (int)Math.round((double)nCols/2);
							if (this.alignToStart)
							{
								colsStart = midPoint;
							}
							else if (this.alignToEnd)
							{
								colsEnd = midPoint;
							}

						}

						//sortVals = this.fStatMats[i][ind].sortRowsExponential(-0.3333, 10, colsStart, colsEnd);
						if ((colsEnd-colsStart) < 20)
						{
							sortVals = this.fStatMats[i][ind].sortRowsExponential(-0.6666, 2, colsStart, colsEnd);
						}
						else
						{
							sortVals = this.fStatMats[i][ind].sortRowsExponential(-2.0, 10, colsStart, colsEnd);
						}
					}
					else
					{
						this.fStatMats[i][ind].sortRowsByList(sortVals);
					}
				}
				

				String tablePrefix = tablePrefixes.get(i);
				if (!this.nohtml)
				{
					if (this.readCounts)
					{
						writer.printf("<H4>ReadCounts %s (%d features)</H4>\n", tablePrefix, fStatMats[i][ind].numFeats());
						writer.println(this.fStatMats[i][0].htmlChart(!this.combineStrands, false, false));
					}
					else
					{
						writer.printf("<H4>%s (%d features)</H4>\n", tablePrefix, fStatMats[i][ind].numFeats());
						writer.println(this.fStatMats[i][2].htmlChart(!this.combineStrands, true, true));
					}
				}

				if (!this.nomatlab)
				{
					if (this.readCounts)
					{
						PrintWriter alignmentWriter = new PrintWriter(new FileOutputStream(
								String.format("%s.%s.featType%d.flank%d.extend%d.readCounts.alignments.csv", outputPrefix, tablePrefixes.get(i), onFeatType, this.flankSize, this.extendRead)));
						this.fStatMats[i][0].matlabCsv(alignmentWriter, !this.combineStrands, this.includeCoords);
						alignmentWriter.close();
					}
					else
					{
						PrintWriter alignmentWriter = new PrintWriter(new FileOutputStream(
								String.format("%s.%s.featType%d.flank%d.extend%d.alignments.csv", outputPrefix, tablePrefixes.get(i), onFeatType, this.flankSize, this.extendRead)));
						this.fStatMats[i][2].matlabCsv(alignmentWriter, !this.combineStrands, this.includeCoords);
						alignmentWriter.close();

						alignmentWriter = new PrintWriter(new FileOutputStream(
								String.format("%s.%s.featType%d.flank%d.extend%d.alignmentWeights.csv", outputPrefix, tablePrefixes.get(i), onFeatType, this.flankSize, this.extendRead)));
						this.fStatMats[i][3].matlabCsv(alignmentWriter, !this.combineStrands, this.includeCoords);
						alignmentWriter.close();
					}
					
					// Coordinates are not valid if sorted
					if (this.nosort)
					{
						PrintWriter coordWriter = new PrintWriter(new FileOutputStream(
								String.format("%s.flank%d.featType%d.coords.csv", outputPrefix, this.flankSize, onFeatType)));
						int coordInd = (this.readCounts) ? 0 : 2;
						this.fStatMats[i][coordInd].coordCsv(coordWriter);
						coordWriter.close();
					}
				}
				
				double[] colorMinMax = {0.0,1.0};
				if (this.heatmap) this.fStatMats[i][ind].launchSwingHeatmap(colorMinMax);
			} // Sample

			if (sortVals != null)
			{
				System.err.println("About to sort expression vals");
				Arrays.sort(sortVals);
				System.err.println("About to print expression vals");
				for (int i = 0; i < sortVals.length; i++)
				{
					sortWriter.println(sortVals[i].toString());
				}
				System.err.println("Done printing expression vals");
			}
			sortWriter.close();
			

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
		System.err.printf("Processing %s, has %d feats\n", chrStr, feats.num_features(chr));
		int featNum = 0;
		FEAT: while (featit.hasNext())
		{
			featNum++;
			
			SimpleGFFRecord rec = (SimpleGFFRecord)featit.next();
			StrandedFeature.Strand featStrand = rec.getStrand();
			String featName = null; // rec.getSeqName();
			int featS = rec.getStart();
			int featE = rec.getEnd();
			
			double sortVal = Double.NaN;
			if (this.sortByExpression != null)
			{
				sortVal = MethylDbUtils.fetchMeanExpression(chrStr, GFFUtils.getGffRecordName(rec), 
						this.sortByExpression);
				if (Double.isNaN(sortVal)) 
				{
					//System.err.println("Feature has no expression value.");
					continue FEAT;
				}
				if (sortVal < 7) 
				{
					System.err.println("Sortval = " + sortVal);
					System.exit(1);
				}
			}

			if (skipUnoriented)
			{
				// Don't use those without orientation
				if (featStrand == StrandedFeature.UNKNOWN) continue;
			}
			else
			{
				if (featStrand == StrandedFeature.UNKNOWN) rec.setStrand(StrandedFeature.POSITIVE);
				featStrand = rec.getStrand();
			}

			GenomicRangeWithRefpoint flankRange = FeatAligner.getAlignmentpointAndFlank(rec, 
					this.flankSize, this.alignToStart, this.alignToEnd, this.censor);
			int alignmentPoint = flankRange.getRefPoint();
			int flankStart = flankRange.getStart() - this.extendRead;
			int flankEnd = flankRange.getEnd() + this.extendRead;
			
			
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(String.format(
					"Fetching coords: censor=%s\talignToStart=%s\tchr=%s\tfeatS=%d\tfeatE=%d\tfeatStrand=%s\talignmentPoint=%d\tflankS=%d\tflankEnd=%d\t\n",
					""+this.censor, ""+this.alignToStart, chrStr, featS, featE, ""+featStrand, alignmentPoint, flankStart, flankEnd));
					

			try 
			{

				// Zero out
				for (int i = 0; i < nS; i++)
				{
					if (this.censor && this.readCounts)
					{
						// Doing every one is way too slow
						this.fStatMats[i][0].zeroOutRange(flankStart, flankEnd, featName, chrStr, alignmentPoint, featStrand, sortVal);
					}
				}


				// We can do this based on Feat iterator or CpG iterator
				if (this.featTableCenters || this.featTableCompleteExtents || this.featTableStarts || this.featTableStartsAndEnds)
				{
					double mLevels[] = new double[nS];
					for (int i = 0; i < nS; i++)
					{
						// Get the overlapping feats
						FeatDbQuerier params = new FeatDbQuerier();
						params.addRangeFilter(chrStr, flankStart,flankEnd);
						params.addFeatFilter(tablePrefixes.get(i));
						FeatIterator targets = new FeatIterator(params);
						while (targets.hasNext())
						{
							GFFRecord targetRec = targets.next();
							//System.err.println("\tFound target rec: " + GFFUtils.gffBetterString(targetRec));

							StrandedFeature.Strand targetStrand = targetRec.getStrand();
							StrandedFeature.Strand alignmentStrand = targetStrand;
							int tS = targetRec.getStart();
							int tE = targetRec.getEnd();
							boolean done = false;
							boolean doStartAndEnd = this.featTableStartsAndEnds || (this.featTableStarts && (targetStrand==StrandedFeature.UNKNOWN));
							for (int p = tS; !done && (p<=tE); p++)
							{
								
								int chromPos = p;
								if (this.featTableCenters)
								{
									chromPos = (int)((tS+tE)/2); // Center
									done = true;
								}
								else if (doStartAndEnd)
								{
									// if we're on p==tS, do start.  If we're on p==(tS+1), do end.  Otherwise we're done
									if (p==tS)
									{
										chromPos = (targetStrand==StrandedFeature.NEGATIVE) ? tE : tS;
										alignmentStrand = (targetStrand==StrandedFeature.UNKNOWN) ? StrandedFeature.POSITIVE : targetStrand;
									}
									else if (p==(tS+1))
									{
										chromPos = (targetStrand==StrandedFeature.NEGATIVE) ? tS : tE;
										alignmentStrand = (targetStrand==StrandedFeature.UNKNOWN) ? StrandedFeature.NEGATIVE : targetStrand.flip();
										done = true;
									}
								}
								else if (this.featTableStarts)
								{
									// If we get here, we know it's oriented
									chromPos = (targetStrand==StrandedFeature.NEGATIVE) ? tE : tS;
									done = true;
								}

								int count = 1;

								this.fStatMats[i][0].addAlignmentPos(
										chromPos,
										(alignmentStrand == StrandedFeature.NEGATIVE) ? 0.0 : count,
												(alignmentStrand == StrandedFeature.NEGATIVE) ? count : 0.0,
														featName, chrStr, alignmentPoint, featStrand, sortVal);
							}
						}
					}
				}
				else
				{

					// Meth
					MethylDbQuerier params = new MethylDbQuerier();
					for (String featFilter : this.featFilters)
					{
						params.addFeatFilter(featFilter,flankSize);
					}
					params.setMinCTreads(this.minCTreads);
					params.setMaxOppstrandAfrac(this.maxOppStrandAfrac);
					params.addRangeFilter(chrStr,flankStart,flankEnd);

					CpgIteratorMultisample cpgit = new CpgIteratorMultisample(params, tablePrefixes);



					while (cpgit.hasNext()) 
					{
						Cpg[] cpgs = cpgit.next();

						int chromPos = cpgs[0].chromPos;
						StrandedFeature.Strand cpgStrand = cpgs[0].getStrand();




						if (this.extendRead > 0)
						{
							int strandMult = 0;
							if (cpgStrand == StrandedFeature.POSITIVE)
							{
								strandMult = 1;
							}
							else if (cpgStrand == StrandedFeature.NEGATIVE)
							{
								strandMult = -1;
							}
							else
							{
								System.err.println("Why did we get  astrandless CpG?");
							}
							chromPos += (  strandMult * this.extendRead );
						}

						// First individual stats
						double mLevels[] = new double[nS];
						for (int i = 0; i < nS; i++)
						{
							double mLevel = cpgs[i].fracMeth(params.getUseNonconversionFilter());
							mLevels[i] = mLevel;
							if (!this.nometh)
							{
								// Meth val
								this.fStatMats[i][2].addAlignmentPos(
										chromPos,
										(cpgStrand == StrandedFeature.NEGATIVE) ? Double.NaN : mLevel,
												(cpgStrand == StrandedFeature.NEGATIVE) ? mLevel: Double.NaN,
														featName, chrStr, alignmentPoint, featStrand, sortVal);
								// And the weights
								this.fStatMats[i][3].addAlignmentPos(
										chromPos,
										(cpgStrand == StrandedFeature.NEGATIVE) ? Double.NaN : cpgs[i].getCpgWeight(),
												(cpgStrand == StrandedFeature.NEGATIVE) ? cpgs[i].getCpgWeight() : Double.NaN,
														featName, chrStr, alignmentPoint, featStrand, sortVal);
							}

							double count = (readCountsIgnoreTotalReads) ? 1.0 : cpgs[i].totalReads;
							if (this.readCounts)
							{
								this.fStatMats[i][0].addAlignmentPos(
										chromPos,
										(cpgStrand == StrandedFeature.NEGATIVE) ? 0.0 : count,
												(cpgStrand == StrandedFeature.NEGATIVE) ? count : 0.0,
														featName, chrStr, alignmentPoint, featStrand, sortVal);
							}
						}

					}


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
