package edu.usc.epigenome.genomeLibs;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.*;

public class ChromScoresBasepairIntersections extends ChromScoresArrayInt {

	
	private Map<String,Integer> f_idmap = new HashMap<String,Integer>();
	private Map<Integer,String> f_digitmap = new HashMap<Integer,String>();
	
	public ChromScoresBasepairIntersections(String genome) {
		super(genome);
	}

	protected int[] getVector(String chr, int pos)
	{
		// Should we auto this to 0?
		int[] ar = (int[])f_chrom_arrays.get(chr);
		return getVectorByBinary(ar[pos]);
	}
	
	
	public int getBinaryById(String id)
	{
		int digit = getDigitById(id);
		return getBinaryByDigit(digit);
	}
	
	
	public int getDigitById(String id)
	{
		Integer digit_obj = (Integer)f_idmap.get(id);
		
		int digit;
		if (digit_obj == null)
		{
			digit = f_idmap.size() + 1;
			f_idmap.put(id, new Integer(digit));
			f_digitmap.put(new Integer(digit),id);
		}
		else
		{
			digit = digit_obj.intValue();
		}
		
		return digit;
	}
	
	public String getIdByDigit(int digit)
	{
		return (String)f_digitmap.get(new Integer(digit));
	}
	
	public String[] Ids()
	{
		String [] ids = new String[f_idmap.size()];
		for (int digit=0; digit < f_idmap.size(); digit++)
		{
				ids[digit] = getIdByDigit(digit+1);
		}
		return ids;
	}
	
	public String IdStr()
	{
		String []ids =  Ids();
		return ListUtils.excelLine(ids);
	}
	
	/********** Binary logic *********/
	
	static protected int getBinaryByDigit(int digit)
	{
		return (int)Math.pow(2,digit-1);
	}
	
	
	protected int[] getVectorByBinary(int binary)
	{
		int[] out = new int[f_idmap.size()];

		for (int digit=1; digit <= f_idmap.size(); digit++)
		{
			int current_binary = getBinaryByDigit(digit);
			out[digit-1] = ((current_binary & binary) > 0)?1:0;
		}
		
		return out;
	}
	
	
	/******* OUTPUT ********/

	
	public static void outputCountMap(Map<Integer,Integer> map, PrintWriter out)
	{
        // One iteration to get the max		
		Iterator<Integer> it = map.keySet().iterator();
		int mx = 0;
		while (it.hasNext())
		{
			Integer binary = (Integer)it.next();
			mx = Math.max(mx, binary.intValue());
		}
		System.err.println("Mx=" + mx);

		
		it = map.keySet().iterator();
		while (it.hasNext())
		{
			Integer binary = (Integer)it.next();
			Integer count = (Integer)map.get(binary);

			System.err.println("Binary=" + binary + "\tcount=" + count);

			
			int pos = 0;
			int current_binary;
			while ((current_binary = (int)Math.pow(2, pos)) < mx)
			{
				int val = ((current_binary & binary.intValue()) > 0) ? count.intValue() : 0;
				if (pos>0) out.print(",");
				out.print(val);

				pos++;
			}
			
			// Make the final column an "other" row
			out.print(",");
			out.print( (binary.intValue() == 0) ? count.intValue() : 0);

			// And line return
			out.println();
		}
	}
	
	public void addSelfToCountMap(String chr, Map<Integer,Integer> map)
	{
		int[] array = (int[])f_chrom_arrays.get(chr);

		int chr_start = minPos(array);
		int chr_end = maxPos(array);
		
		long bases_explored = 0;
		for (int i = chr_start; i<=chr_end; i++)
		{
			if ((bases_explored++ % 1000000) ==0)
			{
				System.err.println(chr + "\t" + bases_explored + " counts added to map");
			}

			int binary = array[i];
			
			Integer key = new Integer(binary);
			Integer count = (Integer)map.get(key);
			
			if (count == null) count = new Integer(0);
			
			count = new Integer(count.intValue()+1);
			map.put(key, count);
		}
	}	
	
	public void outputArrayChr(String chr, PrintWriter out)
	{
		
		int[] array = (int[])f_chrom_arrays.get(chr);

		int chr_start = minPos(array);
		int chr_end = maxPos(array);
		
		
		long bases_explored = 0;
		for (int i = chr_start; i<=chr_end; i++)
		{
			if ((bases_explored++ % 1000000) ==0)
			{
				System.err.println(chr + "\t" + bases_explored + " array lines written");
			}

			int binary = array[i];
			
			if (binary>0)
			{
				int[] vec = getVectorByBinary(binary);
				
				out.println(ListUtils.excelLine(vec));
			}
		}
	}	
	
	 /****** PARSING  ******
	  * 
	  */
	 public int addFeats(String cse_fn, String target_chr_str)
	 throws Exception
	 {
		 BufferedReader in = new BufferedReader(new FileReader(cse_fn));

		 int on_line = 0;
		 String line;
		 
		 
		 String id = cse_fn;
		 int digit = getDigitById(id);
		 
		 int feats_added = 0;
		 while ((line=in.readLine()) != null)
		 {
			 if ((++on_line % 100000) == 0)
			 {
				 System.err.println("\t(chr" + target_chr_str + ")\tOn line\t" + on_line + "\t" + cse_fn);
			 }

			 String[] line_items = line.split(",");
			 

			 
			 String chr = null;
			 try
			 {
				 chr = "chr" + Integer.parseInt(line_items[0]);
			 }
			 catch (Exception e)
			 {
				System.err.println("Skipping bad chromosome: " + line_items[0]); 
			 }
			 
			 if ((line_items.length==3) && (chr == target_chr_str))
			 {
				 int s = Integer.parseInt(line_items[1]);
				 int e = Integer.parseInt(line_items[2]);
				 
				 addRange(target_chr_str,s,e,new Integer(digit));
				 feats_added++;
			 }
		 }

		 in.close();
		 
		 return feats_added;
	 }
	 

	 
	
	/**** Functions from superclass. ****/

	 protected Object addScoreToArray(Object array, int pos, Number score)
	 {
		 int[] array_int = (int[])array;
		 
		 int digit = ((Integer)score).intValue();
		 int bin = getBinaryByDigit(digit);
		 
	//	 System.err.println("oldval=" + array_int[pos] + "\tadding=" + digit + "\tbinary=" + bin + "\tnew=" + (array_int[pos] | bin));
		 try
		 {
			 array_int[pos] = array_int[pos] | bin;
		 }
		 catch (Exception e)
		 {
			 System.err.println("Illegal position: " + pos + "\n" + e.toString());
		 }

		 return array;
	 }

	 
}



