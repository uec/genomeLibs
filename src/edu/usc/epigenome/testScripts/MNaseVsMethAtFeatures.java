package edu.usc.epigenome.testScripts;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintWriter;
import java.sql.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

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
import edu.usc.epigenome.scripts.BisulfiteConvertFasta;


public class MNaseVsMethAtFeatures {

	private static final String C_USAGE = "Use: MNaseVsMethAtFeatures -flankSize 2000 feats.gtf mnaseDbPrefix methDbPrefix";
	private static final double C_NAN = Double.NaN;
	
    @Option(name="-flankSize",usage="bp flanking each side of the feature center")
    protected int flankSize = 2000;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();


	// class vars
	double[][][] fMnaseMat;  // 2 x nFeats x (2*flankSize+1)
	double[][][] fMethMat; // 3 x nFeats x (2*flankSize+1)
	int fCurFeatInd = 0;
	Connection fConn;


	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception	
	{
		new MNaseVsMethAtFeatures().doMain(args);
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

			if(arguments.size() != 3) {
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
		String mnaseDbPrefix = arguments.get(1);
		String methDbPrefix = arguments.get(2);
	
		// Start the actual work.  Parse GFF file
		ChromFeatures cf = new ChromFeatures(featsFn, true);
		cf = cf.centered_regions(1, null);
		int nCf = cf.num_features(11); // ***** CHANGE THIS IF GOING ACROSS CHROMS *** 
		
		// Create arrays
		int nC = (flankSize*2) + 1;
		fMnaseMat = new double[2][nCf][nC];
		MatUtils.initMat(fMnaseMat, 0.0);
		fMethMat = new double[3][nCf][nC];
		MatUtils.initMat(fMethMat[0], C_NAN);
		MatUtils.initMat(fMethMat[1], C_NAN);
		MatUtils.initMat(fMethMat[2], 0.0);
		
		// Setup DB
		setupDb();

		
		for (int chr = 11; chr <= 11; chr++)
		{
			processChrom(chr, cf, mnaseDbPrefix, methDbPrefix);
		}

		// Finish DB
		cleanupDb();
		
		// Output matrices
		PrintWriter writer;
		
		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.fw.csv", featsFn, mnaseDbPrefix)));
		MatUtils.matlabCsv(writer, this.fMnaseMat[0], fCurFeatInd, 0);
		writer.close();

		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.rev.csv", featsFn, mnaseDbPrefix)));
		MatUtils.matlabCsv(writer, this.fMnaseMat[1], fCurFeatInd, 0);
		writer.close();

		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.nmCpgs.csv", featsFn, methDbPrefix)));
		MatUtils.matlabCsv(writer, this.fMethMat[0], fCurFeatInd, 0);
		writer.close();

		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.nReads.csv", featsFn, methDbPrefix)));
		MatUtils.matlabCsv(writer, this.fMethMat[1], fCurFeatInd, 0);
		writer.close();

		writer = new PrintWriter(new FileOutputStream(String.format("%s.%s.nCpgs.csv", featsFn, methDbPrefix)));
		MatUtils.matlabCsv(writer, this.fMethMat[2], fCurFeatInd, 0);
		writer.close();
}
	
	
	protected void processChrom(int chr, ChromFeatures cf, String mnaseDbPrefix, String methDbPrefix)
	throws Exception
	{
		// Setup the DB queries
		String chrStr = cf.public_chrom_str(chr);

		String sql = String.format("select chromPos, strand FROM %s_%s WHERE chromPos >= ? AND chromPos <= ?;",
				mnaseDbPrefix, chrStr);
		PreparedStatement mnasePrep = fConn.prepareStatement(sql);

		sql = String.format("select chromPos, strand, methylReads, totalReads FROM %s_%s WHERE chromPos >= ? AND chromPos <= ?;",
				methDbPrefix, chrStr);
		System.err.println(sql);
		PreparedStatement methPrep = fConn.prepareStatement(sql);

		
		Iterator cfit = cf.featureIterator(chr);
		System.err.println("Processing " + chrStr);
		while (cfit.hasNext())
		{
			GFFRecord rec = (GFFRecord)cfit.next();
			StrandedFeature.Strand featStrand = rec.getStrand();
			
			// Don't use those without orientation
			if (featStrand == StrandedFeature.UNKNOWN) continue;
				
			//System.err.println("On record: " + GFFUtils.gffCsvLine(rec));
			
			int start = rec.getStart() - flankSize;
			int end = rec.getStart() + flankSize;
			//System.err.printf("region len=%d\tarray size=%d\n", end-start+1, fMnaseMat[0].length);
			
			// Mnase
			mnasePrep.setInt(1, start);
			mnasePrep.setInt(2, end);
			ResultSet rs = mnasePrep.executeQuery();
			while (rs.next()) 
			{
				int chromPos = rs.getInt(1);
				String strand = rs.getString(2);
				int strandInd = (strand.charAt(0) == '-') ? 1 : 0;
				
				// Our strand is *relative* to the feature
				int strandIndRel = (featStrand == StrandedFeature.NEGATIVE) ? (1-strandInd) : strandInd;
				
				int relPos = relPos(start, end, chromPos, featStrand);
				//System.err.printf("\tGot mnase row: %d, %d, %s\n", chromPos, relPos, strand, strandInd);
				
				this.fMnaseMat[strandIndRel][this.fCurFeatInd][relPos] = 1.0; // Even if we have multiple, we count it as 1
			}
			
			// Meth
			methPrep.setInt(1, start);
			methPrep.setInt(2, end);
			rs = methPrep.executeQuery();
			while (rs.next()) 
			{
				int chromPos = rs.getInt(1);
				String strand = rs.getString(2);
				//int strandInd = (strand.charAt(0) == '-') ? 1 : 0;
				int methylReads = rs.getInt(3);
				int reads = rs.getInt(4);
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
	
	protected void setupDb()
	throws Exception
	{
		String connStr = "jdbc:mysql://localhost/nucs?user=benb";
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		System.err.println("Getting connection for " + connStr);
		fConn = DriverManager.getConnection(connStr);
		
	}
	
	protected void cleanupDb()
	throws Exception
	{
		fConn.close();
	}

//	static GFFEntrySet parseGff(String fn)
//	throws Exception
//	{
//		GFFEntrySet out = new GFFEntrySet();
//		GFFParser parser = new GFFParser();
//		parser.parse(new BufferedReader(new FileReader(fn)), out.getAddHandler());
//		
//		return out;
//	}
	
}
