package edu.usc.epigenome.scripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintWriter;
import java.sql.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biojava.bio.program.gff.GFFEntrySet;
import org.biojava.bio.program.gff.GFFParser;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.seq.StrandedFeature;


import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.GFFUtils;
import edu.usc.epigenome.genomeLibs.LocUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgIterator;


public class MethylDbToFeatAlignments {

	private static final String C_USAGE = "Use: MethylDbToFeatAlignments -maxFeatSize 10 -skipUnoriented -flankSize 2000 -outputPrefix outputTag feats.gtf";
	private static final double C_NAN = Double.NaN;
	
	@Option(name="-skipUnoriented",usage="If set, skip any unoriented feature (default false)")
	protected boolean skipUnoriented = false;
    @Option(name="-maxFeatSize",usage="maximum size of features to include (default Inf)")
    protected int maxFeatSize = Integer.MAX_VALUE;
    @Option(name="-flankSize",usage="bp flanking each side of the feature center (default 2000)")
    protected int flankSize = 2000;
    @Option(name="-outputPrefix",usage="Prefix for output files (default methylDb)")
    protected String outputPrefix = "methylDb";
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();


	// class vars
	double[][][] fMethMat; // 3 (methylReads,reads,cpgs) x nFeats x (2*flankSize+1)
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
		int nC = (flankSize*2) + 1;
		fMethMat = new double[3][nFeats][nC];
		MatUtils.initMat(fMethMat[0], C_NAN);
		MatUtils.initMat(fMethMat[1], C_NAN);
		MatUtils.initMat(fMethMat[2], 0.0);
		
		for (int chr = 11; chr <= 11; chr++)
		{
			processChrom(chr, feats, skipUnoriented);
		}

		
		// Output matrices
		PrintWriter writer;

		String featsFnBase = (new File(featsFn)).getName();

		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.nmCpgs.csv", outputPrefix, featsFnBase)));
		MatUtils.matlabCsv(writer, this.fMethMat[0], fCurFeatInd, 0);
		writer.close();

		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.nReads.csv", outputPrefix, featsFnBase)));
		MatUtils.matlabCsv(writer, this.fMethMat[1], fCurFeatInd, 0);
		writer.close();

		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.nCpgs.csv", outputPrefix, featsFnBase)));
		MatUtils.matlabCsv(writer, this.fMethMat[2], fCurFeatInd, 0);
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
			//System.err.printf("region len=%d\tarray size=%d\n", end-start+1, fMnaseMat[0].length);
			
	
			
			// Meth
			CpgIterator cpgit = new CpgIterator(chrStr, start, end); 

			while (cpgit.hasNext()) 
			{
				Cpg cpg = cpgit.next();
				
				int chromPos = cpg.chromPos;
				//String strand = cpg.getStrand();
				int methylReads = cpg.totalReadsC(true);
				int reads = cpg.totalReadsCorT(true);
				int relPos = relPos(start, end, chromPos, featStrand);
				//System.err.printf("\tGot meth row: %d, %d, %s = %d/%d\n", chromPos, relPos, strand, methylReads, reads);
				
				this.fMethMat[0][this.fCurFeatInd][relPos] = (double)methylReads;
				this.fMethMat[1][this.fCurFeatInd][relPos] = (double)reads;
				this.fMethMat[2][this.fCurFeatInd][relPos] = 1.0;
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
