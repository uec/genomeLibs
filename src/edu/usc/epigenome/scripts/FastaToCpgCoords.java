package edu.usc.epigenome.scripts;

import java.io.*;
import java.util.*;

import org.biojava.bio.seq.io.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.io.*;
import org.biojava.bio.seq.*;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

public class FastaToCpgCoords {

	private static String C_USAGE = "Use: FastaToCpgCoords file1.fa file2.fa ...";

//	@Option(name="-cpgToY",usage="if true, CpH are set to T and CpG are set to Y (default false)")
//	private boolean cpgToY = false;
//	@Option(name="-repeatMaskToN",usage="if true, lower case bases are set to n (default false)")
//	private boolean repeatMaskToN = false;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();


	public static void main(String[] args)
	throws Exception
	{
		new FastaToCpgCoords().doMain(args);
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



			if(arguments.size() < 1) {
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

//		Alphabet maskeddna = SoftMaskedAlphabet.getInstance(DNATools.getDNA());
//		SymbolTokenization dnaParser = maskeddna.getTokenization("token");

		for (String fn : arguments)
		{


			// Now read the file
			BufferedReader f_read = new BufferedReader(new FileReader(fn));
			SequenceIterator seqs = SeqIOTools.readFastaDNA(f_read);

			
			// Make output file
			String suffix = ".cpgCoords.csv";
			String out_fn = fn + suffix;
			System.err.println("out=" + out_fn);
			PrintWriter out = new PrintWriter(new FileOutputStream(out_fn));

			long total_len = 0;
			int total_seqs = 0;
			while (seqs.hasNext())
			{
				Sequence seq = seqs.nextSequence();
				String seqName = seq.getName();
				seqName = seqName.replaceAll("chr", "");
				
				// System.err.println("Seq: " + seq.seqString());
				int seq_len = seq.length();
				total_len += seq_len;
				System.err.println("Processing " +  seqName + " total len = " + total_len);


					SymbolList symbols = seq;
					String full_residues = symbols.seqString().toUpperCase(); // They come out as all lower case
					

					for (int i = 0; i < seq_len; i++)
					{
						char c = full_residues.charAt(i);
						char c2 = (i==(seq_len-1)) ? 'N' : full_residues.charAt(i+1);
							
							if ((c == 'C') && (c2 == 'G'))
							{
								// CpG 
								out.println(seqName + "," + (i)); // The C
								out.println(seqName + "," + (i+1)); // The reverse strand C
							}
					}
			}

			System.out.println("seqs: " + total_seqs + "\ttotal length: " + total_len);

			f_read.close();
			out.close();
		}

	}
}


