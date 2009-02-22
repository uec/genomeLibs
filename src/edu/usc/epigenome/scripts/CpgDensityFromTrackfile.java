package edu.usc.epigenome.scripts;

import java.util.Iterator;
import java.util.List;
import java.util.ArrayList;
import java.io.*;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.symbol.SimpleSymbolListFactory;
import org.biojava.bio.symbol.Symbol;
import org.biojavax.bio.seq.RichSequence.IOTools;
import org.biojavax.bio.seq.io.FastaFormat;
import org.biojavax.bio.seq.io.SimpleRichSequenceBuilderFactory;
import org.biojavax.bio.seq.RichSequenceIterator;

import org.kohsuke.args4j.*;
import org.kohsuke.args4j.spi.*;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.TrackFiles.*;

public class CpgDensityFromTrackfile {

	private static final String USAGE = "Usage: CpgDensityFromTrackfile -infile outfile.trackl -chrom chr21 -s 150 -e 200";


	@Option(name="-infile",usage="A binary trackfile formatted file to read")
	private String infile = null;
	@Option(name="-chrom",usage="Chromosome")
	private String chrom = null;
	@Option(name="-s",usage="start coord (1-based)")
	private int s = 0;
	@Option(name="-e",usage="end coord (1-based)")
	private int e = 0;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();

	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception
	{
		new CpgDensityFromTrackfile().doMain(args);
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

			if( infile==null || s==0 || e==0 || chrom==null)
			{
				throw new CmdLineException("Must all parameters");
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
		TrackFile trackFile = new TrackFileRandomAccess(new File(infile), "hg18", false, false);
		
		double mean = trackFile.getValAvg(chrom, s, e-s+1);
		System.err.println("mean\t" + mean);
		
	}




}

