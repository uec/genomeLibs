package edu.usc.epigenome.genomeLibs;

import java.util.*;

import org.biojava.bio.symbol.*;
import org.biojava.bio.program.gff.*;

public class InEntrySetFilter implements GFFRecordFilter {
	
	Set f_set = null;
	boolean f_negative = false;
	
	public InEntrySetFilter(GFFEntrySet other, boolean negative)
	{ 
		f_negative = negative;
		
		f_set = new HashSet();
		Iterator it = other.lineIterator();
		while (it.hasNext())
		{
			f_set.add((GFFRecord)it.next());
		}
	}
	
	
	public boolean accept(GFFRecord record) 
	{
		boolean in_set = f_set.contains(record);

		return (f_negative) ? !in_set : in_set;
	}
}
