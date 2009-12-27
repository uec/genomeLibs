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

import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAligner;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerAveraging;
import edu.usc.epigenome.genomeLibs.FeatAligners.FeatAlignerEachfeat;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;


public class MethylDbToFeatAlignments {

	private static final String C_USAGE = "Use: MethylDbToFeatAlignments -maxFeatSize 10 -skipUnoriented -flankSize 2000 -outputPrefix outputTag feats1.gtf feats2.gtf ...";
	
	@Option(name="-skipUnoriented",usage="If set, skip any unoriented feature (default false)")
	protected boolean skipUnoriented = false;
    @Option(name="-maxFeatSize",usage="maximum size of features to include (default Inf)")
    protected int maxFeatSize = Integer.MAX_VALUE;
    @Option(name="-flankSize",usage="bp flanking each side of the feature center (default 2000)")
    protected int flankSize = 2000;
    @Option(name="-outputPrefix",usage="Prefix for output files (default methylDb)")
    protected String outputPrefix = "methylDb";
    @Option(name="-tablePrefix",usage="Prefix for DB table (default " + CpgIterator.DEFAULT_TABLE_PREFIX + ")")
    protected String tablePrefix = null;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();


	// class vars
//	double[][][] fMethMat; // 3 (methylReads,reads,cpgs) x nFeats x (2*flankSize+1)
	FeatAligner[] fMats = new FeatAligner[4]; // 3 (methylReads,reads,cpgs, mLevel)
	int fCurFeatInd = 0;


	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception	
	{
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).setLevel(Level.SEVERE);
		new MethylDbToFeatAlignments().doMain(args);
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

			if(arguments.size() != 1) {
				System.err.println(C_USAGE);
				System.exit(1);
			}
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

		String featsFn = arguments.get(0);

	
		// Start the actual work.  Parse GFF file
		ChromFeatures feats = new ChromFeatures(featsFn, true);
		feats = feats.filterBySize(0, this.maxFeatSize);
		feats = feats.centered_regions(1, null);
		int nFeats = feats.num_features(11); // ***** CHANGE THIS IF GOING ACROSS CHROMS *** 
		
		// Create arrays
//		this.fMats[0] = new FeatAlignerEachfeat(flankSize, false, nFeats);
//		this.fMats[1] = new FeatAlignerEachfeat(flankSize, false, nFeats);
		this.fMats[0] = new FeatAlignerAveraging(flankSize, false);
		this.fMats[1] = new FeatAlignerAveraging(flankSize, false);
		this.fMats[2] = new FeatAlignerAveraging(flankSize, true);
		this.fMats[3] = new FeatAlignerAveraging(flankSize, false);
	
		for (int chr = 11; chr <= 11; chr++)
		{
			processChrom(chr, feats, skipUnoriented);
		}

		
		// Output matrices
		PrintWriter writer;

		String featsFnBase = (new File(featsFn)).getName();

//		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.nmCpgs.csv", outputPrefix, featsFnBase)));
//		this.fMats[0].matlabCsv(writer, true);
//		writer.close();
//
//		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.nReads.csv", outputPrefix, featsFnBase)));
//		this.fMats[1].matlabCsv(writer, true);
//		writer.close();
//
//		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.nCpgs.csv", outputPrefix, featsFnBase)));
//		this.fMats[2].toAverageFeatAligner().matlabCsv(writer, true);
//		writer.close();
//
//		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.mLevel.csv", outputPrefix, featsFnBase)));
//		this.fMats[3].toAverageFeatAligner().matlabCsv(writer, true);
//		writer.close();

		
		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.charts.html", outputPrefix, featsFnBase)));
		writer.printf("<H1>Sample=%s, Feature=%s</H1>\n<P></P>\n", outputPrefix, featsFnBase);
		for (int i = 1; i >=0; i--)
		{
			boolean strandSpec = (i == 0) ? true : false;

			writer.printf("<H3>mLevel</H3>");
			writer.print(this.fMats[3].htmlChart(strandSpec, true, true));

			writer.printf("<H3>Read counts</H3>");
			writer.print(this.fMats[1].htmlChart(strandSpec, true, false));

//			writer.printf("<H3>mCpGs</H3>");
//			writer.print(this.fMats[0].htmlChart(strandSpec, false, false));

			writer.printf("<H3>CpGs</H3>");
			writer.print(this.fMats[2].htmlChart(strandSpec, false, false));
}
		writer.close();
	
	}
	
	
	protected void processChrom(int chr, ChromFeatures feats, boolean skipUnoriented)
	throws Exception
	{
		// Setup the DB queries
		String chrStr = feats.public_chrom_str(chr);

		Iterator featit = feats.featureIterator(chr);
		System.err.println("Processing " + chrStr);
		while (featit.hasNext())
		{
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
				
			//System.err.println("On record: " + GFFUtils.gffCsvLine(rec));
			
			int start = rec.getStart() - flankSize;
			int end = rec.getStart() + flankSize;
//			//System.err.printf("region len=%d\tarray size=%d\n", end-start+1, fMnaseMat[0].length);

	
			
			// Meth
			CpgIterator cpgit = new CpgIterator(chrStr, start, end, tablePrefix); 

			while (cpgit.hasNext()) 
			{
				Cpg cpg = cpgit.next();
				
				int chromPos = cpg.chromPos;
				StrandedFeature.Strand cpgStrand = cpg.getStrand();
				int methylReads = cpg.totalReadsC(true);
				int reads = cpg.totalReadsCorT(true);
				//int relPos = relPos(start, end, chromPos, featStrand);
				Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).info(
						String.format("\tGot meth row: %d, %s = %d/%d\n", chromPos, cpgStrand, methylReads, reads));

				this.fMats[0].addAlignmentPos(chromPos, 
						(cpgStrand == StrandedFeature.NEGATIVE) ? Double.NaN : (double)methylReads,
								(cpgStrand == StrandedFeature.NEGATIVE) ? (double)methylReads : Double.NaN,
										featName, chrStr, rec.getStart(), featStrand);	
				
				this.fMats[1].addAlignmentPos(chromPos, 
						(cpgStrand == StrandedFeature.NEGATIVE) ? Double.NaN : (double)reads,
								(cpgStrand == StrandedFeature.NEGATIVE) ? (double)reads : Double.NaN,
										featName, chrStr, rec.getStart(), featStrand);	

				this.fMats[2].addAlignmentPos(chromPos, 
						(cpgStrand == StrandedFeature.NEGATIVE) ? Double.NaN : (double)1.0,
								(cpgStrand == StrandedFeature.NEGATIVE) ? (double)1.0 : Double.NaN,
										featName, chrStr, rec.getStart(), featStrand);	

				double mLevel = (double)methylReads/(double)reads;
				//System.err.printf("mLevel = %.2f\n",mLevel*100.0);
				this.fMats[3].addAlignmentPos(chromPos, 
						(cpgStrand == StrandedFeature.NEGATIVE) ? Double.NaN : mLevel,
								(cpgStrand == StrandedFeature.NEGATIVE) ? mLevel: Double.NaN,
										featName, chrStr, rec.getStart(), featStrand);	

			}
			
			// Increment feat ind
			this.fCurFeatInd++;
		}
		
		

	}
	
	public static int relPos(int rangeStart, int rangeEnd, int inPos, StrandedFeature.Strand strand)
	{
		int out;
		
		if (strand == StrandedFeature.NEGATIVE)
		{
			out = rangeEnd - inPos;
		}
		else
		{
			out = inPos - rangeStart;
		}
		
		return out;
	}
	
	
}
