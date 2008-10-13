package edu.usc.epigenome.genomeLibs;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.LinkedList;
import java.text.*;




// Keep coverage in a location object for fast 
//iteration

// Factor out type-specific code so we can
// subclass it with different types


import java.util.*;



abstract public class ChromScoresFast {

	/* Class vars */
//	private final static boolean C_BED_STYLE = true;
	
	/* Obj vars */
	protected Map<String,Object> f_chrom_arrays = new HashMap<String,Object>();
//	private Map f_chrom_locs = new HashMap();
	int f_first_chrom = 0;
	int f_last_chrom = 0;
	String f_genome = null;
	
	/* counters */
	protected int f_ranges_added = 0;
	
	/* Constructors */
	public ChromScoresFast() 
	{
		System.err.println("Initializing ChromScoresFast");
		init();
	}
	
	public ChromScoresFast(String genome) 
	{
		init(genome);
	}

	/* Initialization */
	protected void init(String genome)
	{
		System.err.println("Initializing with genome " + genome);
		f_genome = genome;
	}
	
	protected void init()
	{
		init(null);
	}
	
	
	// Override this to change the underlying data structur
	abstract protected Object addScoreToArray(Object array, int pos, Number score); // Does not "sum" score, it replaces it
	abstract protected Number getArrayScore(Object array, int pos);
	//abstract protected double[] getArrayScores(Object array, int st, int end);
	abstract protected double[] getAllScores(Object array);
	abstract protected int minPos(Object array);
	abstract protected int maxPos(Object array);
	
	
	public String[] activeChroms()
	{
		return (String[])(f_chrom_arrays.keySet().toArray(new String[1]));
	}
	
	public boolean chromosomeEmpty(String chr)
	{
		Object out = f_chrom_arrays.get(chr);
		return (out == null);
	}
	
	private Object checkChrom(String chr)
	{
		Object out = f_chrom_arrays.get(chr);
		if (out == null)
		{
			try
			{
				initChrom(chr);
				out = f_chrom_arrays.get(chr);
			}
			catch (Exception e)
			{
				System.err.println("Could not init chrom " + chr + ": " + e.toString());
			}
		}
//		System.err.println("Checking chrom " + chr + " out=" + out);
		return out;
	}
	
	private void initChrom(String chr)
	throws Exception
	{
		System.err.println("Initializing chrom " + chr + " with genome " + f_genome);
		
		// Get the length of the chrom
		int chr_len = GoldAssembly.chromLengthStatic(chr,f_genome);
		
		// Initialize the chrom array
		Object chrom_arr = this.newChromArray(chr_len);
		f_chrom_arrays.put(chr, chrom_arr);
	}
	
	// Probably override this depending on type
	abstract protected Object newChromArray(int chr_len);
//	{
//		return new float[chr_len];
//	}
	
	/* Managing scores */
	
	
	public void addRange(String chr, int s, int e, Number score)
	{
		//System.err.println("Adding range " + chr + ":" + s + "-" + e + "\tcount=" + count);
		for (int i = s; i <= e; i++)
		{
			addScore(chr, i, score);
		}
		
		f_ranges_added++;
		if ((f_ranges_added % 1000) == 0)
		{
			System.err.println(f_ranges_added + "\tranges added");
		}
				
	}
	
	// String chr is optional
	public void addScore(String chr, int pos, Number score)
	{

		// Check if the chromosome is initialized
		Object chrom_array = this.checkChrom(chr);

		
		// Now add score
		this.addScoreToArray(chrom_array, pos, score);
	}
	


	// String chr is optional
	public Number getScore(String chr, int pos)
	{

		// Check if the chromosome is initialized
		Object chrom_array = this.checkChrom(chr);
		return getArrayScore(chrom_array,pos);
	}

	public double[] getScores(String chr, int st, int end, boolean rev)
	{
		int size = Math.abs(end-st+1);
		double[] out = new double[size];
		
		//System.err.println("St: " + st + " end: " + end);
		
		for (int i = 0; i < size; i++)
		{
			int pos = (rev) ? (end-i) : (i+st);
			double score = ((Number)getScore(chr, pos)).doubleValue();
			//System.err.println("\tscore=" + score);
			out[i]=score;
		}
		
		return out;
	}
	
	/*
	 * Static populators
	 */
	
	
	public void populateFromWig(String fn)
	throws Exception
	{
		BufferedReader stream = new BufferedReader(new FileReader(fn));

		String chr = null;
		int step = 0, pos = 0;
		int step_type = 0; // 1 = fixed, 2 = variable
		
		String line = null;
		int line_num = 0;
		while ((line = stream.readLine()) != null)
		{
			
			if (line.startsWith("fixedStep"))
			{
				String[] tab_secs = line.split("\t");
				for (int i = 0; i<tab_secs.length; i++)
				{
					step_type = 1;
					String[] key_val = tab_secs[i].split("=");
					if (key_val[0].equalsIgnoreCase("chrom"))
					{
						chr = key_val[1];
					}
					else if (key_val[0].equalsIgnoreCase("start"))
					{
						pos = Integer.parseInt(key_val[1]);
					}
					else if (key_val[0].equalsIgnoreCase("step"))
					{
						step = Integer.parseInt(key_val[1]);
					}
				}
			}
			else if (line.startsWith("variableStep"))
			{
				step_type = 2;
			}
			else if (line.startsWith("#") || line.startsWith("track"))
			{
				// Just skip
			}
			else
			{
				if (step_type == 1) // Fixed step
				{
					Double val = Double.valueOf(line);
		//			System.err.println("Adding " + chr + ": " + pos + " = " + val);
					this.addScore(chr, pos, val);
					pos += step;
				}
				else if (step_type == 2) // Variable step
				{
					throw new Exception(fn + ", variableStep type not yet supported");
				}
				else if (step_type == 0)
				{
					throw new Exception(fn + ", No fixedStep or variableStep line preceeding: " + line);
				}
				else
				{
					throw new Exception(fn + ", unknown stepType: " + step_type);
				}
			}

			if ( ((line_num++ % 1E7) == 0) )
			{
				System.err.println("On line " + line_num + " " + fn);
			}

		
		
		}
	}
	
	
	/*
	 * Filters
	 * 
	 */
	
	// Sums the score at each coordinate
	public void addScores(ChromScoresFast other)
	{
		String[] active_chroms = this.activeChroms();
		
		for (int i = 0; i < active_chroms.length; i++)
		{
			String chr = active_chroms[i];
			
			Object this_array = f_chrom_arrays.get(chr);
			Object other_array = other.f_chrom_arrays.get(chr);
			
			int chr_start = Math.min(minPos(this_array),minPos(other_array));
			int chr_end = Math.max(maxPos(this_array),maxPos(other_array));

			
			for (int pos = chr_start; pos <= chr_end; pos++)
			{
				if ((pos % 1000000) == 0)
				{
					System.err.println("Summing window " + pos );
				}
				this.addScoreToArray(this_array, pos, other.getArrayScore(other_array, pos));
			}
		}
	}
	
	// If smooth_window is 0, we just return the original
	public ChromScoresFast smooth(int smooth_window, int read_length)
	throws Exception
	{
		return smooth(smooth_window,Double.MIN_VALUE, read_length);
	}
	
	// If smooth_window is 0, we just return the original
	public ChromScoresFast smooth(int smooth_window, double min_avg, int read_length)
	throws Exception
	{
		if (smooth_window <= 0) return this;
		
		// Make the output the same type as the input
		ChromScoresFast out = new ChromScoresArray();
		String[] active_chroms = this.activeChroms();
		
		for (int i = 0; i < active_chroms.length; i++)
		{
			String chr = active_chroms[i];
			Object chrom_array = f_chrom_arrays.get(chr);
						
			if (chrom_array != null)
			{
				addSmoothedChr(smooth_window, chr, chrom_array, out, min_avg, read_length);
			}
		}
		
		return out;
	}
	
//	abstract protected void addSmoothedChr(int smooth_window, String chr, Object chrom_array, ChromScoresFast out);

	protected void addSmoothedChr(int smooth_window, String chr, Object array, ChromScoresFast out, double min_avg, int read_length)
	{
		
		int smooth_half = Math.round((float)smooth_window/(float)2.0)-1;
		int chr_start = minPos(array); // f_first_coords[chr_num];
		int chr_end = maxPos(array); // f_last_coords[chr_num];
		
		LinkedList<Number> queue = new LinkedList<Number>();
		double queue_size = 0.0;
		double running_total = 0.0;
		for (int wind_start = chr_start; wind_start <= (chr_end-smooth_window+1); wind_start++)
		{
			int wind_end = wind_start + smooth_window - 1;
			int wind_center = wind_start + smooth_half;

			Number next_score_obj = this.getArrayScore(array, wind_end);
			if (next_score_obj == null) next_score_obj = new Double(0.0);
			double next_score_double = next_score_obj.doubleValue();


			if ((wind_start % 1000000) == 0)
			{
				System.err.println("Smoothing window " + wind_start + "/" + chr_end + " (minavg=" + min_avg + ")");
			}
			
			
			// Add the wind end onto the end of the queue
			queue.addLast(next_score_obj);
			
			// Add this to the running total
			running_total += next_score_double;
			queue_size++;
			
			// Remove the beginning of the queue if we're up to capacity
			if (queue_size > smooth_window)
			{
				Number remove_score_obj = (Number)queue.removeFirst();
				double remove_score_double = remove_score_obj.doubleValue();
				running_total -= remove_score_double;
				queue_size--;
				//System.err.println("\t\tRemoving");
			}
		
			// Now take the average. If it's 0, leave it out
			double avg = running_total;
			
			if (read_length > 0)
			{
				avg = avg / (double)read_length;
			}
			else
			{
				avg = avg / (double)queue_size;
			}
			
			if ((avg != 0.0) && (avg>min_avg))
			{
				out.addScore(chr, wind_center, new Double(avg));
			}

//			System.err.println("chrom=" + chr + "\twind_start=" + wind_start + "\twind_end=" + 
//					wind_end + "\tqueue size = " + queue.size() + "avg:" + avg);
			

		}		
	}	

	
	/* 
	 * Output
	 */


	// If compact is false, we put each nucleotide on it's own line.  Otherwise we consolidate adjacent
	// nucleotides that have the same value.
	//
	// FORMAT= 0 (wig full), 1 (wig compact), 2 (GADA), 3 (BED)
	public void countOutput(String name, double min_count, PrintWriter out)
	throws Exception
	{
		Iterator<String> it = GoldAssembly.chromIterator(f_genome, false);
		while (it.hasNext())
		{
			String chr = it.next();

			Object chrom_array = f_chrom_arrays.get(chr);
			
			if (chrom_array != null)
			{
				int s = minPos(chrom_array);
				int e = maxPos(chrom_array);
				System.err.println("Writing " + chr + ":" + s + "-" + e);
				countOutputChr(name, min_count, chr, chrom_array,out);
			}
		}
	}
	
	protected void countOutputChr(String name, double min_score, String chr, Object array, PrintWriter out)
	{
		int chr_start = minPos(array); // f_first_coords[chr_num];
		int chr_end = maxPos(array); // f_last_coords[chr_num];
				
		long bases_explored = 0;
		for (int i = chr_start; i<=chr_end; i++)
		{
			if ((bases_explored++ % 1000000) ==0)
			{
				System.err.println(chr + "\t" + bases_explored + " bases explored");
			}

			Number count = this.getArrayScore(array, i);
			double count_double = 0;
			if (count!=null)
			{
				count_double = count.doubleValue();
			}
			
			if (count_double >= min_score)
			{
				String outstr = count.toString();
				out.println(outstr);
			}
			
		}
	}	
	
	// If compact is false, we put each nucleotide on it's own line.  Otherwise we consolidate adjacent
	// nucleotides that have the same value.
	//
	// FORMAT= 0 (wig full), 1 (wig compact), 2 (GADA), 3 (BED)
	public void wigOutput(String name, String fn, boolean print_header, int format, int step)
	throws Exception
	{
		// What's bigger than NEGATIVE_INFINITY?
		wigOutput(name, -1000000000.0,fn, print_header, format, step);
	}
	
	public void wigOutput(String name, double min_count, String fn, boolean print_header, int format, int step)
	throws Exception
	{
		// Open a file writer
		PrintWriter out = new PrintWriter(new FileOutputStream(fn));
		
		wigOutput(name,min_count, out, print_header, format, step);
		
		// Close file writer 
		out.close();
	}
	
	// If compact is false, we put each nucleotide on it's own line.  Otherwise we consolidate adjacent
	// nucleotides that have the same value.
	//
	// FORMAT= 0 (wig full), 1 (wig compact), 2 (GADA), 3 (BED)
	public void wigOutput(String name, PrintWriter out, boolean print_header, int format, int step)
	throws Exception
	{
		// What's bigger than NEGATIVE_INFINITY?
		wigOutput(name, -1000000000.0, out, print_header, format, step);
	}

	
	// If compact is false, we put each nucleotide on it's own line.  Otherwise we consolidate adjacent
	// nucleotides that have the same value.
	//
	// FORMAT= 0 (wig full), 1 (wig compact), 2 (GADA), 3 (BED)
	public void wigOutput(String name, double min_count, PrintWriter out, boolean print_header, int format, int step)
	throws Exception
	{
		wigOutput(name,min_count,out,print_header,format,step,false, null);
	}
	
	public void wigOutput(String name, double min_count, PrintWriter out, boolean print_header, int format, int step, boolean upside_down, String color)
	throws Exception
	{

		int strand_start = 0; // -1
		int strand_end = 0; // 1
		
		for (int s=strand_start; s <=strand_end ; s++)
		{
			String desc = name + " generated by " + this.getClass();
			String head = "track type=wiggle_0 name=\"" + name + "\" description=\"" + desc +
			"\" visibility=2 maxHeightPixels=128:25:11 graphType=bar windowingFunction=mean autoScale=off";
			
			if (color != null)
			{
				head += " color=" + color;
			}

			// System.err.println("Min value = " + min_count);

			if (print_header) out.println(head);

			Iterator<String> it = GoldAssembly.chromIterator(f_genome, false);
			while (it.hasNext())
			{
				String chr = it.next();

				Object chrom_array = f_chrom_arrays.get(chr);

				if (chrom_array != null)
				{
					int start = minPos(chrom_array);
					int end = maxPos(chrom_array);
					System.err.println("Writing " + chr + ":" + start + "-" + end);
					wigOutputChr(name, min_count, chr, chrom_array,out, format,s, step,  upside_down);
				}
			}
		}
	}


	// Strand must be [-1,0,1] where -1 is minus strand, +1 is plus strand, 0 is both
	//
	// FORMAT= 0 (wig full), 1 (wig compact), 2 (GADA)
	public void wigOutputChr(String name, double min_score, String chr, PrintWriter out, int format, int strand, int step, boolean upside_down)
	{
		Object chrom_array = f_chrom_arrays.get(chr);
		wigOutputChr(name, min_score, chr, chrom_array, out, format, strand, step, upside_down);
	}
	
	protected void wigOutputChr(String name, double min_score, String chr, Object array, PrintWriter out, int format, int strand, boolean upside_down)
	{
		wigOutputChr( name,  min_score,  chr,  array,  out,  format,  strand, 1, upside_down);
	}

	protected void wigOutputChr(String name, double min_score, String chr, Object array, PrintWriter out, int format, int strand, int step, boolean upside_down)
	{
		if (format == 1)
		{
			wigOutputChrCompact(name, min_score, chr, array, out, strand, step,  upside_down);
		}
		else if (format == 0)
		{
			wigOutputChrFull(name, min_score, chr, array, out, strand, step,  upside_down);
		}
		else if (format == 2)
		{
			gadaOutputChr(min_score, chr, array, out, strand,  upside_down);
		}
	}	
	

	// Strand must be [-1,0,1] where -1 is minus strand, +1 is plus strand, 0 is both
	//
	protected void wigOutputChrFull(String name, double min_score, String chr, Object array, PrintWriter out, int strand, boolean upside_down)
	{
		wigOutputChrFull( name,  min_score,  chr,  array,  out,  strand, 1, upside_down);
	}

	protected void wigOutputChrFull(String name, double min_score, String chr, Object array, PrintWriter out, int strand, int step, boolean upside_down)
	{

		int chr_start = minPos(array); // f_first_coords[chr_num];
		int chr_end = maxPos(array); // f_last_coords[chr_num];
		
		out.println("variableStep\tchrom=" + chr);
		
		long bases_explored = 0;
		for (int i = chr_start; i<=chr_end; i+=step)
		{
			if ((bases_explored++ % 1000000) ==0)
			{
				System.err.println(chr + "\t" + bases_explored + " bases explored");
			}

			Number count = this.getArrayScore(array, i);
			double count_double = 0;
			if (count!=null)
			{
				count_double = count.doubleValue();
			}
			
			if (count_double >= min_score)
			{
				if (upside_down) count_double = count_double * -1;
				String outstr = i + "\t" + count_double;
				out.println(outstr);
			}
			
		}
	}
	
	
	
	
	// Strand must be [-1,0,1] where -1 is minus strand, +1 is plus strand, 0 is both
	//
//	abstract protected void wigOutputChr(String name, double min_count, String chr, Object chrom_array, PrintWriter out);
	protected void wigOutputChrCompact(String name, double min_count, String chr, Object array, PrintWriter out, int strand, int step, boolean upside_down)
	{

		int chr_start = minPos(array);
		int chr_end = maxPos(array);
		
		out.println("fixedStep\tchrom=" + chr + "\tstart=" + chr_start + "\tstep=" + step +
				"\tspan=500");
		
		NumberFormat format = new DecimalFormat("0.00");
		
		long bases_explored = 0;
		for (int i = chr_start; i<=chr_end; i+=step)
		{
			if ((bases_explored++ % 1000000) ==0)
			{
				System.err.println(chr + "\t" + bases_explored + " bases explored");
			}

			Number count = this.getArrayScore(array, i);
			double count_double = 0.0;
			if (count!=null)
			{
				count_double = count.doubleValue();
			}
			
//			String outstr = chr + "\t" + i + "\t" + (i+1) + "\t" + count_double;
			String outstr = format.format(count_double); //Double.toString(count_double);
			out.println(outstr);
			
		}
	}	
	
		
	
	// Strand must be [-1,0,1] where -1 is minus strand, +1 is plus strand, 0 is both
	//
	protected void gadaOutputChr(double min_score, String chr, Object array, PrintWriter out, int strand, boolean upside_down)
	{

		int chr_start = minPos(array); // f_first_coords[chr_num];
		int chr_end = maxPos(array); // f_last_coords[chr_num];
		
		if (array.getClass().isAssignableFrom(edu.usc.epigenome.genomeLibs.ChromScoresMap.class))
		{
			
		}
		
		long bases_explored = 0;
		for (int i = chr_start; i<=chr_end; i++)
		{
			if ((bases_explored++ % 1000000) ==0)
			{
				System.err.println(chr + "\t" + bases_explored + " bases explored");
			}

			Number count = this.getArrayScore(array, i);
			double count_double = 0;
			if (count!=null)
			{
				count_double = count.doubleValue();
			}
			
			if ((count_double >= min_score) && (count_double != 0.0))
			{
				long rounded = Math.round(count_double);
				if ((double)rounded == count_double)
				{
					out.println(rounded);
				}
				else
				{
					out.println(count_double);
				}
			}
			
		}
	}	

//	// Averaging genomic features
//	// We use a general double[][] array to store values
//	// ar[i][j] where i is the feature number and j is the position
//	// n_i = pre_gran + core_gran + post_gran
//	//
//	// gran = granularity
//	static public double[][] initFeatureAvsArray(FeatureCountParams params, int n_feats)
//	{
//		int total_out_size = params.pre_gran + params.core_gran + params.post_gran;
//		double[][] out = new double[n_feats][total_out_size];
//		return out;
//	}
//	
//	public double[][] addFeatureCounts(ChromFeatures cfs, FeatureCountParams params)
//	throws Exception
//	{
//		int n_feats = (int)cfs.num_features_genome();
//		double[][] mat = initFeatureAvsArray(params, n_feats);
//		addFeatureCounts(cfs, params, mat, 0, null);
//		return mat;
//	}
//	
//	public void addFeatureCounts(ChromFeatures cfs, FeatureCountParams params, double[][] mat, 
//			int feat_start_ind, GFFEntrySet feats_ordered)
//	throws Exception
//	{
//		Iterator it = cfs.featureIterator();
//		int on_feat = feat_start_ind;
//		while (it.hasNext())
//		{
//			GFFRecord feat = (GFFRecord)it.next();
//			feats_ordered.add(feat);
//			featureCounts(feat, params, mat[on_feat]);
//			on_feat++;
//		}
//	}
//
//	public void featureCounts(GFFRecord feat, FeatureCountParams params, double[] out)
//	throws Exception
//	{
//		String feat_chr = feat.getSeqName();
//		int feat_s = feat.getStart();
//		int feat_e = feat.getEnd();
//		boolean feat_rev = (feat.getStrand() == StrandedFeature.NEGATIVE);
//		
//		
//		int pre_s = 0, pre_e = params.pre_gran - 1;
//		int core_s = pre_e + 1, core_e = core_s + params.core_gran - 1;
//		int post_s = core_e + 1, post_e = post_s + params.post_gran - 1;
//		
//		if (post_e != (out.length-1)) System.err.println("Array size ("+out.length+") != expected(" + post_e+1 + ")");
//		
//		int in_s = 0, in_e = 0;
//		int out_s = 0, out_e = 0;
//		double[] scores;
//
//		// 5'
//		in_s = feat_s - ((feat_rev) ? params.post_size : params.pre_size);
//		in_e = feat_s - 1;
//		out_s = (feat_rev) ? post_s : pre_s;
//		out_e = (feat_rev) ? post_e : pre_e;
//		scores = this.getScores(feat_chr, 0, in_s, in_e, feat_rev);
//		MatUtils.downscaleArray(out, scores, out_s, out_e);
//		
//		// Core
//		in_s = feat_s;
//		in_e = feat_e;
//		out_s = core_s;
//		out_e = core_e;
//		scores = this.getScores(feat_chr, 0, in_s, in_e, feat_rev);
//		MatUtils.downscaleArray(out, scores, out_s, out_e);
//		
//		// 3'
//		in_s = feat_e + 1;
//		in_e = feat_e + ((feat_rev) ? params.pre_size : params.post_size);
//		out_s = (!feat_rev) ? post_s : pre_s;
//		out_e = (!feat_rev) ? post_e : pre_e;
//		scores = this.getScores(feat_chr, 0, in_s, in_e, feat_rev);
//		MatUtils.downscaleArray(out, scores, out_s, out_e);
//	}	
	
//	
//	
////	CHARTING
// 
//	public void histogramChr(String name, double min_score, int chr, int strand)
//	{
//		ChromFeatures cf = new ChromFeatures();
//		cf.setGenomeVersion(this.f_genome);
//		String chr_str = cf.public_chrom_str(chr);	
//		
//		Object chrom_array = f_chrom_arrays.get(chr_str);
//
//		if (chrom_array != null)
//		{
//			double[] series = getAllScores(chrom_array);
//			
//			
//			ChromScoresHistogram hist = new ChromScoresHistogram(name, series, 30);
//			hist.pack();
//			RefineryUtilities.centerFrameOnScreen(hist);
//			hist.setVisible(true);
//		}
//	}
//	

	
}
