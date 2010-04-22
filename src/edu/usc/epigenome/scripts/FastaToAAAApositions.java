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

public class FastaToAAAApositions {

	private static String C_USAGE = "Use: FastaToAAAApositions -n 4  file1.fa file2.fa ...";

	@Option(name="-nA",usage="number of As in row to qualify (default 4)")
	private int nA = 4;
	@Option(name="-omitChr",usage="if true, just output position and strand (default false)")
	private boolean omitChr = false;
	@Option(name="-allowOverlaps",usage="Allow adjacent AAAA oligomers to overlap (default false)")
	private boolean allowOverlaps = false;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();


	public static void main(String[] args)
	throws Exception
	{
		new FastaToAAAApositions().doMain(args);
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
			String out_fn = String.format("%s.polyA.n%d.coords.csv",fn,this.nA);
			System.err.println("out=" + out_fn);
			PrintWriter out = new PrintWriter(new FileOutputStream(out_fn));

			long total_len = 0;
			int total_seqs = 0;
			while (seqs.hasNext())
			{
				Sequence seq = seqs.nextSequence();
				String seqName = seq.getName();
				//seqName = seqName.replaceAll("chr", "");
				
				// System.err.println("Seq: " + seq.seqString());
				int seq_len = seq.length();
				total_len += seq_len;
				System.err.println("Processing " +  seqName + " total len = " + total_len);


					SymbolList symbols = seq;
					String full_residues = symbols.seqString().toUpperCase(); // They come out as all lower case
					

					for (int i = 0; i <= (seq_len-this.nA); i++)
					{
						// Do FW and REV (we can only pass one)
						boolean found = false;
						for (int strand=0; !found && (strand<=1); strand++)
						{
							char matchChar = (strand==0) ? 'A' : 'T';
							
							boolean fails=false;
							for (int j = 0; !fails && (j < this.nA); j++)
							{
								if (full_residues.charAt(i+j) != matchChar) fails = true;
							}
							
							if (!fails)
							{
								// We got one.  Signal it and advance main counter to our last position
								int polyaPos = (strand==0) ? i : (i+this.nA-1);
								char strandStr = (strand==0) ? '+' : '-';
								String chrStr = (this.omitChr) ? "" : (seqName + ",");
								out.printf("%s%d,%s\n", chrStr, polyaPos, strandStr);
								found = true;
								if (!this.allowOverlaps)
								{
									i += (this.nA-1); // Advance past the poly-A, don't overcount
								}
							}
						}
					}
					
					
			}

			System.out.println("seqs: " + total_seqs + "\ttotal length: " + total_len);

			f_read.close();
			out.close();
		}

	}
}


