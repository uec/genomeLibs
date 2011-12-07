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
import edu.usc.epigenome.genomeLibs.MiscUtils;

public class FastaToOligoPositions {

	private static String C_USAGE = "Use: FastaToAAAApositions -oligo GA -minCopies 1  file1.fa file2.fa ...";

	@Option(name="-minCopies",usage="minimum number of copies to output (default 1)")
	private int minCopies = 1;
	@Option(name="-oligo",usage="Oligo sequence (must be specified)")
	private String oligo = null;
	@Option(name="-omitChr",usage="if true, just output position and strand (default false)")
	private boolean omitChr = false;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();


	public static void main(String[] args)
	throws Exception
	{
		new FastaToOligoPositions().doMain(args);
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

			if(this.oligo == null) {
				System.err.println("Must supply an oligo sequence");
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

		this.oligo = this.oligo.toUpperCase();
//		String oligo_rev = MiscUtils.revCompNucStr(this.oligo);
		
		for (String fn : arguments)
		{


			// Now read the file
			BufferedReader f_read = new BufferedReader(new FileReader(fn));
			SequenceIterator seqs = SeqIOTools.readFastaDNA(f_read);

			
			// Make output file
			String out_fn = String.format("%s.poly%s.minCopies%d.coords.csv",fn,this.oligo,this.minCopies);
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
				total_seqs++;
				//System.err.println("Processing " +  seqName + " total len = " + total_len);


				SymbolList symbols = seq;
				String full_residues = symbols.seqString().toUpperCase(); // They come out as all lower case
				String full_residues_rev = MiscUtils.revCompNucStr(full_residues);
				
				for (int direction = 0 ; direction <= 1; direction++)
				{
					
//					String oligo_or = (direction == 0) ? this.oligo : oligo_rev;
					String residues_or = (direction == 0) ? full_residues : full_residues_rev;
					
					
					int nMatchesCur = 0; // >0 if we're in a stretch of matches
					int matchStart = -1;
					for (int i = 0; i <= (seq_len-this.oligo.length()); i++)
					{
						
						boolean matchThis = true;
						for (int matchi = 0; (matchi<this.oligo.length()) && matchThis; matchi++)
						{
							char residuechar = residues_or.charAt(matchi+i);
							char oligochar = this.oligo.charAt(matchi);
							//System.err.printf("%c == %c, matchThis=%s\n",residuechar, oligochar, matchThis);
							matchThis &= (residuechar == oligochar);
						}
						
						if ((nMatchesCur>0) && !matchThis)
						{
							// Exiting a match section, output
							//System.err.printf("\tExiting match of %d, strand-relative pos %d (minCopies=%d)\n", nMatchesCur,matchStart, this.minCopies);
							if (nMatchesCur >= this.minCopies)
							{
								int pos = (direction==0) ? matchStart : (seq_len-i);
								char strandStr = (direction==0) ? '+' : '-';
								String chrStr = (this.omitChr) ? "" : (seqName + ",");
								out.printf("%s%d,%s,%d\n", chrStr, pos, strandStr,nMatchesCur);
							}
						}
						else if ((nMatchesCur==0) && matchThis)
						{
							//System.err.printf("\tEntering match\n");
							// Starting a match
							matchStart = i;
						}
						
						
						if (matchThis)
						{
							nMatchesCur++;
							i += (this.oligo.length() - 1);
						}
						else
						{
							nMatchesCur = 0;
							matchStart = -1;
						}
					}
					
					// Check if we ended within a match
					if ((nMatchesCur>0) && (nMatchesCur>=this.minCopies))
					{
						// Exiting a match section, output
						int pos = (direction==0) ? matchStart : 0;
						char strandStr = (direction==0) ? '+' : '-';
						String chrStr = (this.omitChr) ? "" : (seqName + ",");
						out.printf("%s%d,%s,%d\n", chrStr, pos, strandStr,nMatchesCur);
					}
					
				}
			}

			System.out.println("seqs: " + total_seqs + "\ttotal length: " + total_len);

			f_read.close();
			out.close();
		}

	}
}


