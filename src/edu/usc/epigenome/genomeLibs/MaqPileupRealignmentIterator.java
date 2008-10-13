package edu.usc.epigenome.genomeLibs;

import java.io.IOException;
import java.util.TreeMap;

public class MaqPileupRealignmentIterator extends RealignmentIterator {

	private int f_total_bases_read = 0;
	
	public MaqPileupRealignmentIterator(String fn, AlignmentPosOptions apos) 
	throws IOException {
		super(fn, apos);
		// TODO Auto-generated constructor stub
	}

	@Override
	protected AlignmentPos nextAlignment()
	throws IOException 
	{
		// TODO Auto-generated method stub
		String line = this.f_open_stream.readLine();
		String[] line_items = line.split("\t");

	
		String line_chr = line_items[0];
		int line_pos = Integer.parseInt(line_items[1]);
		char line_ref = line_items[2].charAt(0);
		int line_count = Integer.parseInt(line_items[3]);
		String base_pileup_string = line_items[4];
//		String base_quals = line_items[5];
		String read_positions = line_items[7];
		
		f_total_bases_read += line_count;
		
		// Make the output object
		AlignmentPos ap = null;
		if (f_options.f_track_positions)
		{
//			ap = new AlignmentPosSnpsPositions(line_ref, line_chr, line_pos);
			ap = new AlignmentPosSnps(line_ref, line_chr, line_pos);
		}
		else if (f_options.f_track_snps)
		{
			ap = new AlignmentPosSnps(line_ref, line_chr, line_pos);
		}
		else
		{
			ap = new AlignmentPos(line_ref, line_chr, line_pos);
			ap.addDepth(maqPileupDepth(base_pileup_string,read_positions));
		}
	
		return ap;
	}
	
	
	
	
	// Takes the list of read positions (12,12,17,18,...) and determines
	// a depth count.  Reads starting at identical positions are counted
	// a maximum of max_identical times.  The input positions are not always in 
	// order in the input file.
	// 
	// Set max_identical to 0 for a faster implementation
	protected int[] maqPileupDepth(String snps_in, String positions)
	{
		// snps string
		char[] in_ar = snps_in.toCharArray();
		int[] out = new int[2];
		int len = in_ar.length;
		
		// positions string
		TreeMap<String,Integer> counts = null;
		String[] pos_strings = null;
		if ((positions != null) && (f_options.f_max_identical > 0))
		{
			counts = new TreeMap<String,Integer>();
			pos_strings = positions.split(",");
		}
		
		for (int i = 1; i < len; i++) // First char is an "@" character
		{
			char c = in_ar[i];

			// Get the strand (0=fwd, 1=rev)
			int strand = -1;
			switch (c)
			{
			case ',':	strand=0; break;
			case 'A':	strand=0; break;
			case 'C':	strand=0; break;
			case 'T':	strand=0; break;
			case 'G':	strand=0; break;
			case 'N':	strand=0; break;
			case '.':	strand=1; break;
			case 'a':	strand=1; break;
			case 'c':	strand=1; break;
			case 't':	strand=1; break;
			case 'g':	strand=1; break;
			case 'n':	strand=1; break;
			default:	System.err.println("Illegal maq char: "+c); break;
			}
			
			if (strand>=0) // Make sure it's a valid position
			{
				boolean add = true;
				if ((pos_strings != null) && f_options.f_max_identical > 0)
				{
					// Check if we've reached the maximum
					String key = strand + "__" + pos_strings[i-1]; // This list has no "@" leading character, index is one behind the SNP list
					int val = (counts.get(key) == null) ? 0 : ((Integer)counts.get(key)).intValue();
					val++; // The current one
					add = (val <= f_options.f_max_identical);
					counts.put(key, new Integer(val));
				}
				
				if (add) out[strand]++;
			}
		}
		
//		System.err.println("[" + out[0] +"," +out[1] + "] = depthFromReadPositions(" + 
//				snps_in + ", " + max_identical + ", " + positions + ")");
		return out;
	}

}
