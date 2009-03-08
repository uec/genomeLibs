package edu.usc.epigenome.testScripts;

import java.io.*;
import java.sql.*;
import java.util.*;

import org.biojava.bio.program.gff.*;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.*;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import org.usckeck.genome.*;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.TrackFiles.TrackFile;
import edu.usc.epigenome.genomeLibs.TrackFiles.TrackFileRandomAccess;

public class GtfToDistTrackfile {


	static final String USAGE = "GtfToDistTrackfile -gtf in.gtf -outfile out.track -genome hg18";

	/**
	 * @param args
	 */

	@Option(name="-gtf",usage="A Gtf formatted file to read")
	private String gtf = null;
	@Option(name="-outfile",usage="A binary track format file (will be overwritten!)")
	private String outfile = null;
	@Option(name="-genome",usage="The genome assembly (hg18, mm9, etc.)")
	private String genome = "hg18";
	@Option(name="-useLong",usage="TBA (default true)")
	private boolean useLong = true;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();

	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception
	{
		new GtfToDistTrackfile().doMain(args);
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

			if( gtf==null || outfile==null)
			{
				throw new CmdLineException("Must supply and input and output file");
			}
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



		// Setup output file
		TrackFile trackFile = new TrackFileRandomAccess(new File(outfile), genome, true, true);		

		// Get features
		ChromFeatures feats = new ChromFeatures(gtf, true);
		
		int[] chroms = {11}; // !!!!! CHANGE BACK
//		for (int chrNum : feats.activeChroms())
		for (int chrNum : chroms)
		{
			String chrStr = feats.public_chrom_str(chrNum);
			int[] chromVals = new int[GoldAssembly.chromLengthStatic(chrStr, genome)]; 
			
			
			// Do a FW sweep then a reverse
			System.err.println("Processing chrom " + chrNum);
			setChromValsOnedir(feats, chrStr, chromVals, 1);

			
			// Write out the chrom
			System.err.println("Writing chrom " + chrStr);
			trackFile.writeChrom(chrStr, chromVals);
			System.err.println("\tDone writing chrom " + chrStr);

		}

		
		
		
	}
	
	
	public void setChromValsOnedir(ChromFeatures feats, String chrStr, int[] chromVals, int passNum)
	throws Exception
	{
		
		int chrNum = feats.chrom_from_public_str(chrStr);
		
		Iterator featInt = feats.featureIterator(chrNum, ((passNum == 0) ? new GffStartCoordComparator() : new GffEndCoordComparator()));

		int currentPos = 0;
		GFFRecord curFeat = null;
		GFFRecord nextFeat = (GFFRecord) featInt.next();
		while (++currentPos <= GoldAssembly.chromLengthStatic(chrStr, genome))
		{
			int currentVal = Integer.MAX_VALUE;
			if ((currentPos % 100000) ==0) System.err.println("On " + chrStr + ", pos " + currentPos);
			
			if (curFeat == null)
			{
				// Haven't gotten to the first gene yet, skip.
			}
			if ((currentPos >= curFeat.getStart()) && (currentPos <= curFeat.getEnd()))
			{
				// in the current feature
				currentVal = 0;
			}
			else if ((nextFeat == null) || (currentPos < nextFeat.getStart()))
			{
				// Haven't gotten to the next feat yet (or there isn't one)
				int dirMultiplier = (curFeat.getStrand() == StrandedFeature.NEGATIVE) ? -1 : 1;
				int dist = dirMultiplier * (currentPos - curFeat.getEnd());
				currentVal = dist;
			}
			else if (currentPos == nextFeat.getStart())
			{
				// We've entered the next gene
				currentVal = 0;
				if (featInt.hasNext())
				{
					curFeat = nextFeat;
					nextFeat = (GFFRecord)featInt.next();
				}
			}
			
			// Set the bit
			chromVals[currentPos-1] = currentVal;
		}
		
	}

}
