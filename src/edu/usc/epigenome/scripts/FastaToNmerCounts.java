package edu.usc.epigenome.scripts;
import java.io.*;
import java.util.*;

import org.biojava.bio.seq.io.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.*;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;




public class FastaToNmerCounts {
	
	private static final String USAGE = "Use: FastaToNmerCounts -listAllPerms -nmer 6 -startPosition 0 -bothStrands [file.fa]\n" + 
	"If file.fa is not included, uses stdin.";
	
    @Option(name="-bothStrands",usage="Reverse complement all sequences (default false)")
    private boolean bothStrands = false;
    @Option(name="-listAllPerms",usage="list all permutations of length nmer on output (default false, list only observed nmers)")
    private boolean listAllPerms = false;
    @Option(name="-nmer", usage="nmer length (default 6)")
    private int nmer = 6;
    @Option(name="-startPosition", usage="If set, we only look at nmers starting at position N of read (value from 1 to readLen-nmer), no revcomp (default none)")
    private int startPosition = -1;
    // receives other command line parameters than options
    @Argument
    private List<String> arguments = new ArrayList<String>();	

	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new FastaToNmerCounts().doMain(args);
    }
    
	public void doMain(String[] args)
	throws Exception
	{


		CmdLineParser parser = new CmdLineParser(this);
        // if you have a wider console, you could increase the value;
        // here 80 is also the default
        parser.setUsageWidth(80);
        BufferedReader f_read;
        try
		{
        	parser.parseArgument(args);


    		if( arguments.isEmpty() )
        	{
        		//throw new CmdLineException("Must supply at least one input file");
    			f_read = new BufferedReader(new InputStreamReader(System.in));
    			System.err.println("Reading from stdin");
        	}
    		else
    		{
        		String fn = arguments.get(0);
        		f_read = new BufferedReader(new FileReader(fn));
  			
    		}

    		if (this.startPosition>=0)
    		{
    			if (this.bothStrands) throw new CmdLineException("Can not set -startPosition and -bothStrands together");
    			if (this.startPosition==0) throw new CmdLineException("-startPosition must be greater than 0");
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




		// Get permutations
		Vector alph = new Vector();
		alph.add("A"); alph.add("C"); alph.add("T"); alph.add("G");
		Vector perms = new Vector();
		getPermutations(perms, alph, nmer, 0, "");
		System.err.println("Found " + perms.size() + " perms");


		// Make counts
		SequenceIterator seqs = SeqIOTools.readFastaDNA(f_read);


		HashMap counts = new HashMap(1 + (4^nmer));
		HashMap onemers = new HashMap(4);
		int on_seq = 1;
		long total_len = 0;
		long total_len_onemers = 0;
		int total_seqs = 0;
		
		int lastPos = (this.startPosition>0) ? (this.startPosition+nmer-1) : 0;
		while (seqs.hasNext())
		{
			Sequence seq = seqs.nextSequence();
			
			addCounts(nmer, counts, seq, this.startPosition);
			if (!this.listAllPerms) addCounts(1, onemers, seq,this.startPosition,lastPos);
			total_len += ((this.startPosition>0) ? 1 : (seq.length()-nmer+1));
			total_len_onemers += ((this.startPosition>0) ? nmer : seq.length());
			if (bothStrands)
			{
				SymbolList revseq = DNATools.reverseComplement(seq);
				addCounts(nmer, counts, revseq, this.startPosition);
				if (!this.listAllPerms) addCounts(1, onemers, revseq,this.startPosition,lastPos);
				total_len += ((this.startPosition>0) ? 1 : (seq.length()-nmer+1));
				total_len_onemers += ((this.startPosition>0) ? nmer : seq.length());
			}
			
			total_seqs++;
			if ( (seq.length() > 1E6) || ((total_seqs % 1E5)==0) )
			{
				System.err.println("Finished seq\t" + total_seqs);
			}
		}
		
		System.out.println("seqs: " + total_seqs + "\ttotal length: " + total_len + "\tonemers=" + total_len_onemers);

		

		if (listAllPerms)
		{
			// Output all permutations
			TreeSet nmers = new TreeSet(perms);
			Iterator it = nmers.iterator();
			while (it.hasNext())
			{
				String nmer = (String)it.next();
				int count;
				if (counts.containsKey(nmer))
				{
					count = ((Integer)counts.get(nmer)).intValue();
				}
				else
				{
					count = 0;
				}
				//count++; // tiny pseudocount to avoid counts of 0
				System.out.println(nmer + "\t" + count);
			}
		}
		else
		{

			// Reverse tree and output permutations in order of 
			// count
			TreeMap revcounts = new TreeMap();
			Iterator it = counts.keySet().iterator();
			while (it.hasNext())
			{
				String residues = (String)it.next();
				Integer count = (Integer)counts.get(residues);
				Integer negcount = new Integer( -1 * count.intValue()); // Neg transformation for descending order

				// We have to make it a set since there may be more than one
				// nmer for a given count
				TreeSet residue_set = (TreeSet)revcounts.get(negcount);
				if (residue_set == null) residue_set = new TreeSet();
				residue_set.add(residues);
				
				revcounts.put(negcount, residue_set);
			}
			it = revcounts.keySet().iterator();
			while (it.hasNext())
			{
				Integer negcount = (Integer)it.next();
				TreeSet residue_set = (TreeSet)revcounts.get(negcount);
				
				Iterator res_it = residue_set.iterator();
				while (res_it.hasNext())
				{
					String residues = (String)res_it.next();
					
					int count = (-1 * negcount.intValue());
					double frac = (double)count / total_len;
					
					double expected = fracExpected(residues,onemers, total_len_onemers);
					double enriched = frac / expected;
					
					System.out.println(residues + "\t" + count + "\t" + frac + "\t" + expected + "\t" + enriched);
					//System.out.println("preAlignmentNmers" + "," + residues + ",+,-1,-1," + count); 
				}
			}
		}


		
		f_read.close();
		
	}
	
	private static double fracExpected(String residues, Map onemers, long total_len)
	{
		
		double expected = 1.0;
		
		for (int i = 0; i < residues.length(); i++)
		{
			String nuc_i = residues.substring(i,i+1);
			//System.err.println(residues + ") Looking for " + nuc_i + " onemer, got " + onemers.get(nuc_i));
			
			Integer count_i = (Integer)onemers.get(nuc_i);
			double expected_i = count_i.doubleValue() / (double)total_len;
			expected *= expected_i;
		}

		return expected;
	}
	
	private static void getPermutations(Vector out_vect, Vector alph, int total_len, int cur_num_done, String cur_str)
	{
		// Done condition
		if (cur_num_done == total_len)
		{
			out_vect.add(cur_str);
		}
		else
		{
			Iterator it = alph.iterator();
			while (it.hasNext())
			{
				String next = (String)it.next();
				getPermutations(out_vect, alph, total_len, cur_num_done+1, cur_str + next);
			}
		}
	}
	
	
	private static void addCounts(int n, Map counts, SymbolList seq, int inStartPosition)
	{
		addCounts(n,counts,seq,inStartPosition,-1);
	}
	
	private static void addCounts(int n, Map counts, SymbolList seq, int inStartPosition, int inEndPosition)
	{
		
		//	int len = seq.length();
		
		String full_residues = seq.seqString();
		//System.err.println("Full str: "+full_residues);
		int len = full_residues.length();
		
		int lastPos = (inEndPosition>0) ? inEndPosition : (len-n+1);
		
		int start = (inStartPosition>0) ? (inStartPosition) : 1;
		int end = ((inStartPosition>=1) && (inEndPosition<1)) ? start : lastPos;
//		System.err.printf("Adding counts\tn=%d,\tstart=%d\tend=%d\n",n,start,end);
		for (int i = start; i <= end; i++)
		{
			int to = i+n-1;
			
			// String residues = seq.subStr(i, to); // This has a huge memory leak
			String residues = full_residues.substring(i-1, to);
			
			
			residues = residues.toUpperCase();
			
			if ((i % 1000000) == 0)
			{
				System.out.println("\tpos=" + i + "-" + to + ": " + residues);
			}
			
			if (residues.indexOf('N') >= 0)
			{
			}
			else if (counts.containsKey(residues))
			{
				Integer last = (Integer)counts.get(residues);
				Integer next = new Integer(1 + last.intValue());
				counts.put(residues, next);
			}
			else
			{
				counts.put(residues, new Integer(1));
			}
		}
	}
  }

