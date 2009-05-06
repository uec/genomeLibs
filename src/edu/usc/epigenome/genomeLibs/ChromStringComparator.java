package edu.usc.epigenome.genomeLibs;


import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ChromStringComparator implements Comparator<String> 
{
	static protected Pattern pat = Pattern.compile("^c(hr)?([0-9]+)$",Pattern.CASE_INSENSITIVE); 

	public int compare(String o1, String o2) {
	
		Matcher m1 = pat.matcher(o1);
		Matcher m2 = pat.matcher(o2);
		
		boolean o1isInt = m1.matches();
		boolean o2isInt = m2.matches();
		
//		System.err.println("Comparing " + o1 + " to " + o2 + "\to1int=" + o1isInt + "\to2int=" + o2isInt);
		
		if (o1isInt && o2isInt)
		{
			Integer o1int = Integer.parseInt(m1.group(2));
			Integer o2int = Integer.parseInt(m2.group(2));
//			System.err.println("\to1int=" + o1int + "\to2int=" + o2int);

			return o1int.compareTo(o2int);
		}
		else if (o1isInt && !o2isInt)
		{
			return -1;
		}
		else if (!o1isInt && o2isInt)
		{
			return 1;
		}
		else //if (!o1isInt && !o2isInt)
		{
			return o1.compareTo(o2);
		}

		
		
	}

	

}
