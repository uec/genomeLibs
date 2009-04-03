package edu.usc.epigenome.scripts;
import java.io.*;
import java.util.*;

import org.biojava.bio.seq.io.*;
import org.biojava.bio.symbol.*;
import org.biojava.bio.seq.io.*;
import org.biojava.bio.seq.*;




public class NmerCounts {
	
	private static String C_USAGE = "Use: NmerCounts n both_strands(0|1) file.fa";
	
	
	static final boolean C_ALL_PERMS = false;
	
	public static void main(String[] args)
	throws Exception {

		if(args.length != 3) {
			System.err.println(C_USAGE);
			System.exit(1);
		}

		int n = Integer.parseInt(args[0]);
		
		int both_strands_in = Integer.parseInt(args[1]);
		boolean both_strands = (both_strands_in>0);

		String fn = args[2];

		// Get permutations
		Vector alph = new Vector();
		alph.add("A"); alph.add("C"); alph.add("T"); alph.add("G");
		Vector perms = new Vector();
		getPermutations(perms, alph, n, 0, "");
		System.err.println("Found " + perms.size() + " perms");


		// Make counts
		BufferedReader f_read = new BufferedReader(new FileReader(fn));
		SequenceIterator seqs = SeqIOTools.readFastaDNA(f_read);


		HashMap counts = new HashMap(1 + (4^n));
		HashMap onemers = new HashMap(4);
		int on_seq = 1;
		long total_len = 0;
		int total_seqs = 0;
		while (seqs.hasNext())
		{
			Sequence seq = seqs.nextSequence();
			
			addCounts(n, counts, seq);
			addCounts(1, onemers, seq);
			total_len += seq.length();
			if (both_strands)
			{
				SymbolList revseq = DNATools.reverseComplement(seq);
				addCounts(n, counts, revseq);
				addCounts(1, onemers, revseq);
				total_len += seq.length();
			}
			
			total_seqs++;
			if ( (seq.length() > 1E6) || ((on_seq % 1E5)==0) )
			{
				System.err.println("Finished seq\t" + on_seq++);
			}
		}
		
		System.out.println("seqs: " + total_seqs + "\ttotal length: " + total_len);

		

		if (C_ALL_PERMS)
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
				count++; // tiny pseudocount to avoid counts of 0
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
					
					double expected = fracExpected(residues,onemers, total_len);
					double enriched = frac / expected;
					
					System.err.println(residues + "\t" + count + "\t" + frac + "\t" + expected + "\t" + enriched);
					System.out.println("preAlignmentNmers" + "," + residues + ",+,-1,-1," + count); 
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
	
	
	private static void addCounts(int n, Map counts, SymbolList seq)
	{
		
		//	int len = seq.length();
		
		String full_residues = seq.seqString();
		//System.err.println("Full str: "+full_residues);
		int len = full_residues.length();
		
		int end = len-n+1;
		for (int i = 1; i <= end; i++)
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

