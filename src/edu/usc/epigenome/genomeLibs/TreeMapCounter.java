/**
 * 
 */
package edu.usc.epigenome.genomeLibs;

import java.util.*;


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
	}
	

	/**
	 * Initialization and termination
	 */

	public void init() {
		// Initialize maps
		//System.err.println("Initializing SymbolCounts");
	}


	public void finish() {
		// Nothing to do
		// System.err.println("Finishing SymbolCounts");
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
	
	public void increment(T key)
	{
		increment(key,1);
	}
	
	public void increment(T key, int numToAdd)
	{
		int newCount = this.getCount(key) + 1;
		this.setCount(key, newCount);
	}
	
	public void setCount(T key, int newCount)
	{
		this.put(key, new Integer(newCount));
	}
	
	
	
	/******  OUTPUT ********/
	
	
	public String excelOutput()
	{
		String out = "";
		
		Iterator<T> strIt = this.keySet().iterator();
		while (strIt.hasNext())
		{
			T key = strIt.next();
			out += key.toString() + "," + this.get(key);
			out += "\n";
		}
		
		return out;
	}

}
