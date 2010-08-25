package edu.usc.epigenome.genomeLibs.ChromScores;

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
import java.util.logging.Logger;

import org.apache.commons.math.linear.MatrixUtils;

import sun.tools.tree.ThisExpression;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.genomeLibs.WigOptions;



abstract public class ChromScoresFast {

	/* Class vars */
//	private final static boolean C_BED_STYLE = true;
	
	/* Obj vars */
	protected Map<String,Object> f_chrom_arrays = new TreeMap<String,Object>();
//	private Map f_chrom_locs = new HashMap();
	int f_first_chrom = 0;
	int f_last_chrom = 0;
	int f_arbitrary_genome_length = (int)100000000;
	String f_genome = null;
	
	public final static String ARBITRARY_GENOME = "arbitraryGenome";
	
	
	/* counters */
	protected int f_ranges_added = 0;
	
	/* Constructors */
//	public ChromScoresFast() 
//	{
//		System.err.println("Initializing ChromScoresFast");
//		init();
//	}
	
	// If genome = ARBITRARY_GENOME, we all arbitrary chromosomes with a set length
	// (settable by arbitraryGenomeLength)
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
	abstract protected Object addScoreToArray(Object array, int pos, Number score); // Increments score
	abstract protected void setScore(Object array, int pos, Number score); // Sets score
	abstract protected Number getArrayScore(Object array, int pos);
	//abstract protected double[] getArrayScores(Object array, int st, int end);
	abstract protected double[] getAllScores(Object array);
	abstract protected int minPos(Object array);
	abstract protected int maxPos(Object array);
	
	
	/**
	 * @return the f_arbitrary_genome_length
	 */
	public int getArbitraryGenomeLength() {
		return f_arbitrary_genome_length;
	}

	/**
	 * @param f_arbitrary_genome_length the f_arbitrary_genome_length to set
	 */
	public void setArbitraryGenomeLength(int in_length) {
		this.f_arbitrary_genome_length = in_length;
	}

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
				System.exit(0);
			}
		}
//		System.err.println("Checking chrom " + chr + " out=" + out);
		return out;
	}
	
	private void initChrom(String chr)
	throws Exception
	{
		if (f_genome == null)
		{
			(new Exception()).printStackTrace();
			throw new Exception("Called ChromScoresFast::initChrom without specifying genome");
		}
		
		// System.err.println("Initializing chrom " + chr + " with genome " + f_genome);
		
		// Get the length of the chrom
		int chr_len = (f_genome == ARBITRARY_GENOME) ? this.getArbitraryGenomeLength() : GoldAssembly.chromLengthStatic(chr,f_genome);
		
		// Initialize the chrom array
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("About to initialize ChromScoresFast");
		Object chrom_arr = this.newChromArray(chr_len);
		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).severe("Done initializing ChromScoresFast");
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
//		int chr_len = 0;
//		try {
//			chr_len = (f_genome == ARBITRARY_GENOME) ? this.getArbitraryGenomeLength() : GoldAssembly.chromLengthStatic(chr,f_genome);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//		if (pos <= chr_len)
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

	public double getScoresTotal(String chr, int st, int end)
	{
		double[] scores = getScores(chr,st,end,false);
		//System.err.printf("Got %d scores\n",scores.length);
		return MatUtils.nanSum(scores);
	}

	public double getScoresTotal(String chr)
	{
		int s = this.chromMinPos(chr);
		int e = this.chromMaxPos(chr);
			
		double total=0.0;
		for (int i = s; i <= e; i++)
		{
			int pos = i;
			double score = ((Number)getScore(chr, pos)).doubleValue();
			//System.err.println("\tscore=" + score);
			total+=score;
		}
		
		return total;
	}
	
	public double getCoverageTotal(String chr)
	{
		int s = this.chromMinPos(chr);
		int e = this.chromMaxPos(chr);
			
		double total=0.0;
		for (int i = s; i <= e; i++)
		{
			int pos = i;
			double score = ((Number)getScore(chr, pos)).doubleValue();
			//System.err.println("\tscore=" + score);
			if (score!=0.0) total+=1.0;
		}
		
		return total;
	}

	


	public int chromMinPos(String chr)
	{
		Object chrom_array = this.checkChrom(chr);
		return this.minPos(chrom_array);
	}
	
	public int chromMaxPos(String chr)
	{
		Object chrom_array = this.checkChrom(chr);
		return this.maxPos(chrom_array);
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
	
	
	/**
	 * @param other A mask array.  Any position set to 0 in the mask array is set to 0 in the primary array.
	 * @param invert If set, any position set to 0 in the mask array is kept
	 * @return
	 */
	public void mask(ChromScoresFast other, boolean invert)
	{
		for (String chr : this.activeChroms())
		{
			int s = this.chromMinPos(chr);
			int e = this.chromMaxPos(chr);
			Object thisChrArr = this.f_chrom_arrays.get(chr);

			for (int i = s; i <= e; i++)
			{
				int pos = i;
				double otherScore = ((Number)other.getScore(chr, pos)).doubleValue();
				boolean toss = (invert) ? (otherScore != 0.0) : (otherScore == 0.0);
				if (toss) this.setScore(thisChrArr, pos, 0.0);
			}
		}
	}
	
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
		return smooth(smooth_window, min_avg, read_length, 1);
	}
	
	public ChromScoresFast smooth(int smooth_window, double min_avg, int read_length, int step)
	throws Exception
	{
		if (smooth_window <= 0) return this;
		
		// Make the output the same type as the input
		ChromScoresFast out = new ChromScoresArray(f_genome);
		out.setArbitraryGenomeLength(this.getArbitraryGenomeLength() / step);
		String[] active_chroms = this.activeChroms();
		
		for (int i = 0; i < active_chroms.length; i++)
		{
			String chr = active_chroms[i];
			Object chrom_array = f_chrom_arrays.get(chr);
						
			if (chrom_array != null)
			{
				addSmoothedChr(smooth_window, chr, chrom_array, out, min_avg, read_length, step);
			}
		}
		
		return out;
	}
	
//	abstract protected void addSmoothedChr(int smooth_window, String chr, Object chrom_array, ChromScoresFast out);

	protected void addSmoothedChr(int smooth_window, String chr, Object array, ChromScoresFast out, double min_avg, int read_length, int step)
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
			
			if ((avg != 0.0) && (avg>min_avg) && ((wind_center%step)==0))
			{
				out.addScore(chr, wind_center/step, new Double(avg));
			}

//			System.err.println("chrom=" + chr + "\twind_start=" + wind_start + "\twind_end=" + 
//					wind_end + "\tqueue size = " + queue.size() + "avg:" + avg);
			

		}		
	}	

	
	/* 
	 * Output
	 */

	public void singleLinePerChrOutput(PrintWriter out)
	throws Exception
	{
		ListUtils.setDelim(",");
		Iterator<String> it = (this.f_genome==ARBITRARY_GENOME) ? this.f_chrom_arrays.keySet().iterator() : GoldAssembly.chromIterator(f_genome, false);
		while (it.hasNext())
		{
			String chr = it.next();
			Object chrom_array = f_chrom_arrays.get(chr);

			double[] scores = (double[])chrom_array; // this.getAllScores(chrom_array);
			// System.err.println("Outputting " + chr + "\tlen=" + scores.length);
			out.println(chr + "," + ListUtils.excelLine(scores));
		}
		out.flush();
	}
	
	// If compact is false, we put each nucleotide on it's own line.  Otherwise we consolidate adjacent
	// nucleotides that have the same value.
	//
	// FORMAT= 0 (wig full), 1 (wig compact), 2 (GADA), 3 (BED)
	public void countOutput(String name, double min_count, PrintWriter out)
	throws Exception
	{
		Iterator<String> it = (this.f_genome==ARBITRARY_GENOME) ? this.f_chrom_arrays.keySet().iterator() : GoldAssembly.chromIterator(f_genome, false);
		
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
	
	public void wigOutput(String fn, WigOptions wo)
	throws Exception
	{
		// Open a file writer
		PrintWriter out = new PrintWriter(new FileOutputStream(fn));
		
		wigOutput(out, wo);
		
		// Close file writer 
		out.close();
	}
	

	
	public void wigOutput(PrintWriter out, WigOptions wo)
	throws Exception
	{

		int strand_start = 0; // -1
		int strand_end = 0; // 1
		
		for (int s=strand_start; s <=strand_end ; s++)
		{
			String head = wo.trackHead();
			if (wo.f_print_head) out.println(head);

			Iterator<String> it = (this.f_genome==ARBITRARY_GENOME) ? this.f_chrom_arrays.keySet().iterator() : GoldAssembly.chromIterator(f_genome, false);
			while (it.hasNext())
			{
				String chr = it.next();

				Object chrom_array = f_chrom_arrays.get(chr);

				if (chrom_array != null)
				{
					int start = minPos(chrom_array);
					int end = maxPos(chrom_array);
					System.err.println("Writing " + chr + ":" + start + "-" + end);
					wigOutputChr(chr, chrom_array, out, wo);
				}
			}
		}
	}


	// Strand must be [-1,0,1] where -1 is minus strand, +1 is plus strand, 0 is both
	//
	// FORMAT= 0 (wig full), 1 (wig compact), 2 (GADA)
	protected void wigOutputChr(String chr, Object array, PrintWriter out, WigOptions wo)
	{
		if (wo.f_format == 1)
		{
			wigOutputChrCompact(chr, array, out, wo);
		}
		else if (wo.f_format == 0)
		{
			wigOutputChrFull(chr, array, out, wo);
		}
		else if (wo.f_format == 2)
		{
			gadaOutputChr(chr, array, out, wo);
		}
	}	
	

	protected void wigOutputChrFull(String chr, Object array, PrintWriter out, WigOptions wo)
	{

		int chr_start = minPos(array); // f_first_coords[chr_num];
		int chr_end = maxPos(array); // f_last_coords[chr_num];
		
		out.println("variableStep\tchrom=" + chr);
		
		long bases_explored = 0;
		for (int i = chr_start; i<=chr_end; i+=wo.f_step)
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
			
			if (count_double >= wo.f_min_score)
			{
				if (wo.f_upside_down) count_double = count_double * -1;
				String outstr = i + "\t" + count_double;
				out.println(outstr);
			}
			
		}
	}
	
	
	
	
	// Strand must be [-1,0,1] where -1 is minus strand, +1 is plus strand, 0 is both
	//
//	abstract protected void wigOutputChr(String name, double min_count, String chr, Object chrom_array, PrintWriter out);
	protected void wigOutputChrCompact(String chr, Object array, PrintWriter out, WigOptions wo)
	{

		int chr_start = minPos(array);
		int chr_end = maxPos(array);
		
		out.print(wo.fixedStepHead(chr,chr_start));
		out.println();
		
		NumberFormat format = new DecimalFormat("0.00");
		
		long bases_explored = 0;
		for (int i = chr_start; i<=chr_end; i+=wo.f_step)
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
			if (wo.f_upside_down) count_double = count_double * -1;
			
			
//			String outstr = chr + "\t" + i + "\t" + (i+1) + "\t" + count_double;
			String outstr = format.format(count_double); //Double.toString(count_double);
			out.println(outstr);
			
		}
	}	
	
		
	
	// Strand must be [-1,0,1] where -1 is minus strand, +1 is plus strand, 0 is both
	//
	protected void gadaOutputChr(String chr, Object array, PrintWriter out, WigOptions wo)
	{

		int chr_start = minPos(array); // f_first_coords[chr_num];
		int chr_end = maxPos(array); // f_last_coords[chr_num];
		
		if (array.getClass().isAssignableFrom(edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresMap.class))
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
			
			if ((count_double >= wo.f_min_score) && (count_double != 0.0))
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


	
}
