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

public class FastaToCpgFile {

	private static final String USAGE = "Usage: FastaToCpgFile -useLong -infile file1.fasta -outfile outfile.trackl";


	@Option(name="-infile",usage="A FASTA formatted file to read")
	private String infile = null;
	@Option(name="-outfile",usage="A binary track format file (will be overwritten!)")
	private String outfile = null;
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
		new FastaToCpgFile().doMain(args);
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

			if( infile==null || outfile==null)
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

		BufferedReader br = null;
		try
		{
			br = new BufferedReader(new FileReader(infile));
		}
		catch (Exception e)
		{
			System.err.println("Could not open input file " + infile);
			System.exit(1);
		}

		// Setup output file
		TrackFile trackFile = new TrackFileRandomAccess(new File(outfile), "hg18", true, true);

		//		// This whole IOTools version seems to use Gb of RAM even after the first few chromosomes
		//		RichSequenceIterator stream = IOTools.readFasta(br, IOTools.getDNAParser(), 
		//                new SimpleRichSequenceBuilderFactory(new SimpleSymbolListFactory()), null); // This one had a memory leak

		// Read input file.
		SequenceIterator stream = SeqIOTools.readFastaDNA(br);


		while (stream.hasNext())
		{
			Sequence seq = (Sequence)stream.nextSequence();
//			String seqStr = seq.seqString(); // Get the seq string, because the seq.symbolAt() function can be very slow


			// Can we find the chrom?
			String chr = seq.getName();
			int expLength = GoldAssembly.chromLengthStatic(chr, "hg18");


			if (expLength != seq.length())
			{
				System.err.println("Error: Chrom " + chr + " exp length (" + expLength +") != actual length (" + seq.length() + ")");
			}
			else
			{
				int[] buf = new int[expLength];

				Symbol lastSymbol = seq.symbolAt(1);
				buf[0] = Integer.MIN_VALUE;
				for (int i = 1; i < expLength; i++)
				{
					Symbol curSymbol = seq.symbolAt(i);
					boolean cpg = ((lastSymbol == DNATools.c()) && (curSymbol == DNATools.g()));
					buf[i] = (cpg) ? 1 : 0;

					if (cpg && ((i%100000)==0)) System.err.println("Found CpG: " + chr + ": " + i);

					lastSymbol = curSymbol;
				}


				System.err.println("Writing " + chr);
				trackFile.writeChrom(chr, buf);
			}
		}
	}





}

