package edu.usc.epigenome.scripts;

import java.io.*;
import java.util.*;

import org.biojava.bio.seq.io.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.io.*;
import org.biojava.bio.seq.*;

public class BisulfiteConvertFasta {

	private static String C_USAGE = "Use: BisulfiteConvert file1.fa file2.fa ...";


	public static void main(String[] args)
	throws Exception {

		if(args.length < 1) {
			System.err.println(C_USAGE);
			System.exit(1);
		}

		for (String fn : args)
		{

			
			// Now read the file
			BufferedReader f_read = new BufferedReader(new FileReader(fn));
			SequenceIterator seqs = SeqIOTools.readFastaDNA(f_read);

			// Make output file 
			String fw_fn = fn;
			fw_fn = fw_fn.replaceFirst(".fa$",".BSfw.fa");
			System.err.println("FW=" + fw_fn);
			PrintWriter fw = new PrintWriter(new FileOutputStream(fw_fn));
			String rev_fn = fn;
			rev_fn = rev_fn.replaceFirst(".fa$",".BSrev.fa");
			System.err.println("REV=" + rev_fn);
			PrintWriter rev = new PrintWriter(new FileOutputStream(rev_fn));

			int on_seq = 1;
			long total_len = 0;
			int total_seqs = 0;
			while (seqs.hasNext())
			{
				Sequence seq = seqs.nextSequence();
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
					String full_residues = symbols.seqString();

					int line_count = 0;
					for (int i = 0; i < seq_len; i++)
					{
						char c = full_residues.charAt(i);

						if (s==0)
						{
							if (c == 'C')
							{
								c = 'T';
							}
							else if (c == 'c')
							{
								c = 't';
							}
						}
						else
						{
							if (c == 'G')
							{
								c = 'A';
							}
							else if (c == 'g')
							{
								c = 'a';
							}
						}

						out.print(c);
						if (line_count++ == 60)
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


