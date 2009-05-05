package edu.usc.epigenome.genomeLibs;


import java.util.Comparator;

public class ChromStringComparator implements Comparator<String> 
{


	public int compare(String o1, String o2) {
	
		o1 = o1.replace("chr", "");
		o2 = o2.replace("chr", "");
		
		
		boolean o1isInt = o1.matches("^[0-9]+$");
		boolean o2isInt = o2.matches("^[0-9]+$");
		
//		System.err.println("Comparing " + o1 + " to " + o2 + "\to1int=" + o1isInt + "\to2int=" + o2isInt);
		
		if (o1isInt && !o2isInt)
		{
			return -1;
		}
		else if (!o1isInt && o2isInt)
		{
			return 1;
		}
		else if (!o1isInt && !o2isInt)
		{
			return o1.compareTo(o2);
		}
		else
		{
			Integer o1int = Integer.parseInt(o1);
			Integer o2int = Integer.parseInt(o2);
			return o1int.compareTo(o2int);
		}
		
		
	}

	

}
