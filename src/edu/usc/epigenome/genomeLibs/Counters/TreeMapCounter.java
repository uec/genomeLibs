/**
 * 
 */
package edu.usc.epigenome.genomeLibs.Counters;

import java.util.*;
import java.io.*;


/**
 * @author benb
 * 
 *
 */
public class TreeMapCounter<T extends Comparable<T>> extends TreeMap<T,Integer>  {

	
	private static final long serialVersionUID = -2401288295505078770L;

	/**
	 * Constructors
	 */
	public TreeMapCounter() 
	{
		super();
	}
	
	
	/**
	 * Accessors
	 */

	public int getCount(T key)
	{
		Integer intObj = this.get(key);
		if (intObj == null) return 0;
		return intObj.intValue();
	}
	
	public int getTotalCount()
	{
		int total = 0;
		for(Integer val : this.values())
		{
			total += val.intValue();
		}
		return total;
	}
	
	public double getFrac(T key)
	{
		double keyCount = (double)getCount(key);
		double total = (double)getTotalCount();
		return (keyCount/total);
	}

	public void increment(T key)
	{
		increment(key,1);
	}
	
	public void decrement(T key)
	{
		decrement(key,1);
	}

	public void increment(Collection<T> keys)
	{
		for(T key : keys)
		{
			this.increment(key);
		}
	}
	
	public void increment(T key, int numToAdd)
	{
			int newCount = this.getCount(key) + numToAdd;
			if (newCount!=0) this.setCount(key, newCount); // If it's still 0, don't create a new key
	}
	
	public void decrement(T key, int numToSubtract)
	{
			int newCount = this.getCount(key) - numToSubtract;
			if (newCount!=0) this.setCount(key, newCount);  // If it's still 0, don't create a new key
	}
	
	public void setCount(T key, int newCount)
	{
		this.put(key, new Integer(newCount));
	}
	
	/**
	 * Adds all key counts from another counter into this one
	 * @param other The tree map to add to this one
	 */
	public void addCounts(TreeMapCounter<T> other)
	{
		for(T tkey : other.keySet())
		{
			this.increment(tkey, other.getCount(tkey));
		}
	}
	
	/**
	 * Adds all key counts from a set of counters into this one
	 * @param maps The tree maps to add to this one
	 */
	public void addCounts(Collection<TreeMapCounter<T>> maps)
	{
		for(TreeMapCounter<T> map : maps)
		{
			addCounts(map);
		}
	}
	
	/* 
	 * This is really a static one but we need it at runtime to determine type of T
	 */
	public TreeMapCounter<T> meanCounter(Collection<TreeMapCounter<T>> counters)
	{
		// Get the totals.
		TreeMapCounter<T> totalCounts = new TreeMapCounter<T>();
		totalCounts.addCounts(counters);
		
		// Divide all by the total number (rounding)
		for (T key : totalCounts.keySet())
		{
			int newVal = Math.round((float)totalCounts.getCount(key) / (float)counters.size());
			totalCounts.setCount(key, newVal);
		}		
		
		return totalCounts;
	}
	
	/******  OUTPUT ********/
	

	public String oneLineOutput()
	{
		String out = "";
		for (T key : this.keySet())
		{
			out += key + "=" + this.getCount(key) + ", ";
		}
		return out;
	}
	
	public String excelOutput()
	{
		return excelOutput(null);
	}

	public String excelOutput(String firstCol)
	{
		OutputStream os = new ByteArrayOutputStream(this.size()*25);
		PrintStream ps = new PrintStream(os);
		this.excelOutput(firstCol, ps);
		return os.toString();
	}
	
	public void excelOutput(String firstCol, PrintStream ps)
	{
		Iterator<T> strIt = this.keySet().iterator();
		while (strIt.hasNext())
		{
			T key = strIt.next();
			if (firstCol!=null) ps.print(firstCol + ",");
			ps.print(key.toString() + "," + this.get(key));
			ps.print("\n");
		}

	}
}
