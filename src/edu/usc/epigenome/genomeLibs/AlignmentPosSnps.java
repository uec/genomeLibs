package edu.usc.epigenome.genomeLibs;

import java.io.*;
import java.util.*;

import org.biojava.bio.program.gff.*;
import org.biojava.bio.seq.StrandedFeature;

public class AlignmentPosSnps extends AlignmentPos {

	/* Class vars */
	public static final double NO_COVERAGE = -1.0;

	/* Obj vars */
	public int[] f_consensus_readpos =  {0,0}; // If -1 , no consensus.  If 0, no reads yet
	public Map f_counts = new TreeMap();// new int[4][2]; // 0=T, 1=C, 2=G, 3=A .   0=fw, 1=rev

	private static Map c_char_map = null;
	private static char[] c_char_map_rev = null;
	private static int[] c_char_map_revcomp_inds = null;
	
	private static Vector c_alph_vect = null;
	private static int c_alph_len = -1;

	static
	{
		initMap();
	}

	/*****************
	 *  Constructors
	 */
	public AlignmentPosSnps(char ref, String chr, int pos)
	{
		initMap();

		f_ref = ref;
		f_chr = chr;
		f_pos = pos;
		
		// Set up a non-position specific matrix
		int[][] mat = new int[4][2]; // 0=A, 1=C, 2=G, 3=T.   0=fw, 1=rev
		Integer key = new Integer(0); // This is the universal one
		f_counts.put(key,mat);
	}

	public AlignmentPosSnps(char ref, String chr, int pos, int[] depth)
	{
		initMap();

		f_ref = ref;
		f_chr = chr;
		f_pos = pos;
		f_depth = depth;
	}
	
	/****************
	 * Utility
	 */

	private static void initMap()
	{
		if (c_char_map == null);
		{	
			c_char_map = new HashMap();
			c_char_map.put("A", new Integer(0));
			c_char_map.put("a", new Integer(0));
			c_char_map.put("C", new Integer(1));
			c_char_map.put("c", new Integer(1));
			c_char_map.put("G", new Integer(2));
			c_char_map.put("g", new Integer(2));
			c_char_map.put("T", new Integer(3));
			c_char_map.put("t", new Integer(3));
		}
		
		if (c_char_map_rev == null)
		{
			c_char_map_rev = new char[] { 'A', 'C', 'T', 'G' };
		}

		if (c_char_map_revcomp_inds == null)
		{
			c_char_map_revcomp_inds = new int[] { 3, 2, 1, 0 };
		}
	
	}

	public static int revCompMapInd(int ind)
	{
		return c_char_map_revcomp_inds[ind];
	}
	
	public static int getMapInd(char c)
	throws Exception
	{
		//System.err.println(c + " " + c_char_map);
		Integer ind = (Integer)c_char_map.get(String.valueOf(c));

		int out = -1;
		if (ind != null)
		{
			out = ind.intValue();
		}
		else if ((c == 'N') || (c == 'n'))
		{
	//		System.err.println("Found N character .. skipping");
		}
		else
		{
			throw new Exception("Can't map value " + c);
		}

		return out;
	}

	public static Iterator alphIterator()
	{
		if (c_alph_vect == null)
		{
			c_alph_vect = new Vector();
			c_alph_vect.add("A");
			c_alph_vect.add("C");
			c_alph_vect.add("G");
			c_alph_vect.add("T");
		}
		return c_alph_vect.iterator();
	}
	
	public static int alphLen()
	{
		if (c_alph_len < 0)
		{
			Iterator it = alphIterator();
			c_alph_len = 0;
			while (it.hasNext())
			{
				it.next();
				c_alph_len++;
			}
		}
		
		return c_alph_len;
	}

	

	/**
	 * 
	 * Getters
	 * 
	 */

	public Iterator readPositions()
	{
		return f_counts.keySet().iterator();
	}


	
	// c should always be the base sequenced relative to the
	// forward strand on the reference sequence.  So if the
	// read is on the opposite strand, you must complement the
	// base.
	
	public int getCount(char c, boolean reference_forward_strand, boolean read_same_strand)
	throws Exception
	{
		return getCount(c, reference_forward_strand,read_same_strand, -1);
	}

	public int getCount(char c, boolean reference_forward_strand, boolean read_same_strand, int read_position)
	throws Exception
	{
		int[][] mat = this.counts(read_position);
		
		if (mat==null)
		{
			return 0;
		}
		else
		{
			boolean read_forward_strand = (reference_forward_strand == read_same_strand);
			if (!reference_forward_strand) c = MiscUtils.revCompNuc(c);
			
			int read_strand_ind = (read_forward_strand) ? 0 : 1;
			int c_ind = AlignmentPosSnps.getMapInd(c);
			return mat[c_ind][read_strand_ind];
		}
	}
	

	
	public int[] getDepth()
	{
		if (f_depth == null)
		{
			f_depth = new int[2];
			f_depth[0] = getDepthFromCounts(true,true);
			f_depth[1] = getDepthFromCounts(true,false);
		}
		
		return f_depth;
	}
	

	
	public int getDepthFromCounts(boolean reference_forward_strand, boolean read_same_strand)
	{
		return getDepth(-1,reference_forward_strand, read_same_strand);
	}
	
	// fw_strand=true returns fw strand depth, false returns rev strand depth
	public int getDepth(int read_pos, boolean reference_forward_strand, boolean read_same_strand)
	{
		int[][] mat = counts(read_pos);

		boolean read_forward_strand = (reference_forward_strand == read_same_strand);
		int read_strand_ind = (read_forward_strand) ? 0 : 1;

		
		int c_len = mat.length;
		int total = 0;
		for (int i = 0; i < c_len; i++)
		{	
			total += mat[i][read_strand_ind];
		}
		return total;
	}	
	
	public double getFrac(char c, boolean reference_forward_strand, boolean read_same_strand)
	throws Exception
	{
		int count =  getCount(c, reference_forward_strand,read_same_strand, -1);
		int depth = getDepth(reference_forward_strand,read_same_strand);
		
		if (count > depth)
		{
			System.err.println("Why is count " + count + " bigger than depth " + depth + "\n" + MatUtils.matString(counts()));
		}

		double frac = NO_COVERAGE;
		if (depth > 0)
		{
			frac = (double)count/(double)depth;
		}
		
		return frac;
	}

	public int[][] counts()
	{
		return counts(-1);
	}
	
	public double[][] fracs(int read_pos)
	{
		double[][] out = new double[4][2]; // 0=A, 1=C, 2=G, 3=T.   0=fw, 1=rev;
		
		int[][] counts = this.counts(read_pos);
		int[] totals = new int[2];
		for (int i = 0; i < 4; i++)
		{
			totals[0] += counts[i][0];
			totals[1] += counts[i][1];
		}
		
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				out[i][j] = (totals[j] == 0) ? 0 : 
					((double)counts[i][j] / (double)totals[j]);
			}
		}		
		
		return out;
	}
	
	public double[][] counts_double(int read_pos)
	{
		double[][] out = new double[4][2]; // 0=A, 1=C, 2=G, 3=T.   0=fw, 1=rev;
		
		int[][] counts = this.counts(read_pos);
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				out[i][j] = (double)counts[i][j];
			}
		}		
		
		return out;
	}

	
	public int[][] counts(int read_pos)
	{
		int[][] out = null; // 0=A, 1=C, 2=G, 3=T.   0=fw, 1=rev;
		
		if (read_pos > 0)
		{
			Integer key = new Integer(read_pos);
		
			// Check if the read position exists
			if (f_counts.containsKey(key))
			{
				out = (int[][])f_counts.get(key);
			}
			else
			{
				out = new int[4][2];
			}
		}
		else
		{
			Iterator it = f_counts.values().iterator();
			while (it.hasNext())
			{

				int [][]mat = (int [][])it.next();
				if (out == null)
				{
					out = mat;
				}
				else
				{
					out = MatUtils.sumMats(out,mat);
				}
				
			}
		}
		
		return out;
	}
	

	
	public AlignmentPos clone(boolean flip_strand)
	{
		char ref = (flip_strand) ? MiscUtils.revCompNuc(this.f_ref) : this.f_ref;
		AlignmentPosSnps ap = new AlignmentPosSnps(ref, this.f_chr, this.f_pos);

		if (f_depth != null)
		{
			ap.f_depth = new int[2];
			ap.f_depth[0] = this.f_depth[1];
			ap.f_depth[1] = this.f_depth[0];
		}

		Map new_counts = new TreeMap();// new int[4][2]; // 0=T, 1=C, 2=G, 3=A .   0=fw, 1=rev
		Iterator it = this.f_counts.keySet().iterator();
		while (it.hasNext())
		{
			Object key = (Object)it.next();
			int [][]mat = (int [][])this.f_counts.get(key);
			//System.err.println("Key=" + key);
			
			// *** BE VERY CAREFUL, THIS INVERSION ONLY WORKS IF ALPH ORDERING
			// IS SET UP THAT WAY ******
			new_counts.put(key, MatUtils.invertedMat(mat));
		}
		ap.f_counts = new_counts;
		
		return ap;
	}
	

	/*****************
	 * 
	 * Setters
	 * 
	 */
	
	
	public void add(char c, boolean read_forward_strand)
	throws Exception
	{
		add(c,read_forward_strand,0);
	}
	
	public void add(char c, boolean read_forward_strand, int read_pos)
	throws Exception
	{
		int strand_ind = (read_forward_strand) ? 0 : 1;
		int map_ind = getMapInd(c);

		if (map_ind >= 0)
		{
			
			// The old way is to make the key the position, if we want to 
			// recover all.
			int read_pos_key = 0; // read_pos 
			Integer key = new Integer(read_pos_key);
			int[][] mat = (int[][])f_counts.get(key);
			if (mat == null)
			{
				mat = new int[4][2]; // 0=A, 1=C, 2=G, 3=T.   0=fw, 1=rev
				f_counts.put(key,mat);
			}



			mat[map_ind][strand_ind]++;

			// Keep track of consensus read positions
			if (f_consensus_readpos[strand_ind] == 0)
			{
				f_consensus_readpos[strand_ind] = read_pos;
			}
			else if (f_consensus_readpos[strand_ind] == -1)
			{
				// Already no consensus
			}
			else 
			{
				if (f_consensus_readpos[strand_ind] == read_pos)
				{
					// do nothing, it's still a consensus
				}
				else
				{
					// No consensus
					f_consensus_readpos[strand_ind] = -1;
				}
			}
//			System.err.println("read_pos=" + read_pos + "consensus=" + f_consensus_readpos);
		}
	}
	


	// Returns the char that was added.
	
	public char addMaqPileupChar(char c)
	throws Exception
	{
		return addMaqPileupChar(c,0);
	}
	
	public char addMaqPileupChar(char c, int read_pos)
	throws Exception
	{
		// Check for maq special characters
		char readstrand_c;
		boolean readstrand_forward_strand;
		
		switch(c)
		{
		case ',': readstrand_c = f_ref; readstrand_forward_strand = true; break;
		case '.': readstrand_c = f_ref; readstrand_forward_strand = false; break;
		default: readstrand_c = c; readstrand_forward_strand = Character.isUpperCase(c); break;
		}


		add(readstrand_c, readstrand_forward_strand, read_pos);
		return readstrand_c;
	}
	
	public void addDepth(int[] depth)
	{
		getDepth(); // Make sure it's set
		f_depth[0] += depth[0];
		f_depth[1] += depth[1];
	}

	
	public void resetCounts()
	{
		this.f_depth = new int[] {0,0};
		this.f_counts = new TreeMap();
	}
	
	



}
