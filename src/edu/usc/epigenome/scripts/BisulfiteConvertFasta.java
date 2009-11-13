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

public class BisulfiteConvertFasta {

	private static String C_USAGE = "Use: BisulfiteConvertFasta -cpgToY -repeatMaskToN file1.fa file2.fa ...";

	@Option(name="-cpgToY",usage="if true, CpH are set to T and CpG are set to Y (default false)")
	private boolean cpgToY = false;
	@Option(name="-repeatMaskToN",usage="if true, lower case bases are set to n (default false)")
	private boolean repeatMaskToN = false;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();


	public static void main(String[] args)
	throws Exception
	{
		new BisulfiteConvertFasta().doMain(args);
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

		Alphabet maskeddna = SoftMaskedAlphabet.getInstance(DNATools.getDNA());
		SymbolTokenization dnaParser = maskeddna.getTokenization("token");

		for (String fn : arguments)
		{


			// Now read the file
			BufferedReader f_read = new BufferedReader(new FileReader(fn));
			RichSequenceIterator seqs = RichSequence.IOTools.readFasta(f_read, dnaParser, null);

			
			// Make output file
			String suffix = (repeatMaskToN) ? ".masked.fa" : ".fa";
			String fw_fn = fn;
			fw_fn = fw_fn.replaceFirst(".fa$",".BSfw" + suffix);
			System.err.println("FW=" + fw_fn);
			PrintWriter fw = new PrintWriter(new FileOutputStream(fw_fn));
			String rev_fn = fn;
			rev_fn = rev_fn.replaceFirst(".fa$",".BSrev" + suffix);
			System.err.println("REV=" + rev_fn);
			PrintWriter rev = new PrintWriter(new FileOutputStream(rev_fn));

			long total_len = 0;
			int total_seqs = 0;
			while (seqs.hasNext())
			{
				Sequence seq = seqs.nextSequence();
				// System.err.println("Seq: " + seq.seqString());
				int seq_len = seq.length();
				total_len += seq_len;
				System.err.println("Processing " + seq.getName() + " total len = " + total_len);

				for (int s = 0; s <= 1; s++)	
				{
					// s==0 means fw
					PrintWriter out = (s==0) ? fw : rev;

					out.print(">" + seq.getName());
					out.println();

					SymbolList symbols = seq;
					String full_residues = symbols.seqString(); // They come out as all lower case
					

					int line_count = 0;
					for (int i = 0; i < seq_len; i++)
					{
						char c = full_residues.charAt(i);

						if (repeatMaskToN && Character.isLowerCase(c))
						{
							c = 'n';
						}
						if (s==0) // s==0 means fw
						{
							char c2 = (i==(seq_len-1)) ? 'N' : full_residues.charAt(i+1);
							
							if (c == 'C')
							{
								c = (cpgToY && (c2 == 'G' || c2 == 'g')) ? 'Y' : 'T'; 
							}
							else if (c == 'c')
							{
								c = (cpgToY && (c2 == 'G' || c2 == 'g')) ? 'y' : 't';
							}
						}
						else
						{
							char c0 = (i==0) ? 'N' : full_residues.charAt(i-1);

							if (c == 'G')
							{
								c = (cpgToY && (c0 == 'C' || c0 == 'c')) ? 'R' : 'A';
							}
							else if (c == 'g')
							{
								c = (cpgToY && (c0 == 'C' || c0 == 'c')) ? 'r' : 'a';
							}
						}

						out.print(c);
						if (line_count++ == 70)
						{
							out.println();
							line_count = 0;
						}
					}

					out.println();
				}
			}

			System.out.println("seqs: " + total_seqs + "\ttotal length: " + total_len);

			f_read.close();
			fw.close();
			rev.close();
		}

	}
}


