package edu.usc.epigenome.genomeLibs;

import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.SimpleGFFRecord;


public class MiscUtils {


	public static void reverseArray(Object[] b) {
		   int left  = 0;          // index of leftmost element
		   int right = b.length-1; // index of rightmost element
		  
		   while (left < right) {
		      // exchange the left and right elements
		      Object temp = b[left]; 
		      b[left]  = b[right]; 
		      b[right] = temp;
		     
		      // move the bounds toward the center
		      left++;
		      right--;
		   }
		}//endmethod reverse
	
	public static String revString(String in)
	{
		StringBuilder sb = new StringBuilder(in);
		sb.reverse();
		return sb.toString();
	}

	
	public static String revCompNucStr(String in)
	{
		String out = "";
		
		for (int i = (in.length()-1); i >= 0; i--)
		{
			out += revCompNuc(in.charAt(i));
		}

		//System.err.println("revCompStr(" + in + ") = " + out);
		return out;
	}

	public static char revCompNuc(char c)
	{
		char cout = c;

		switch (c) {
		case 'A': cout = 'T'; break;
		case 'a': cout = 't'; break;
		case 'T': cout = 'A'; break;
		case 't': cout = 'a'; break;
		case 'C': cout = 'G'; break;
		case 'c': cout = 'g'; break;
		case 'G': cout = 'C'; break;
		case 'g': cout = 'c'; break;

		case 'w': cout = 'w'; break;
		case 'W': cout = 'W'; break;
		case 's': cout = 's'; break;
		case 'S': cout = 'S'; break;
		case 'Y': cout = 'R'; break;
		case 'y': cout = 'r'; break;
		case 'R': cout = 'Y'; break;
		case 'r': cout = 'y'; break;
		}

		return cout;
	}
	
	public static char toggleCase(char c)
	{
		
		char out;
		if (Character.isUpperCase(c))
		{
			out = Character.toLowerCase(c);
		}
		else if (Character.isLowerCase(c))
		{
			out = Character.toUpperCase(c);
		}
		else
		{
			out = c;
		}
		return out;
	}
	
	
	public static String twodArrayString(int[][] a)
	{
		int ni = a.length;
		int nj = (a[0].length);
		
		String out = "";
		for (int i = 0; i < ni; i++)
		{
			for (int j = 0; j < nj; j++)
			{
				if (j>0) out += "\t";
				out += a[i][j];
			}
			out += "\n";
		}
		
		return out;
	}
	
	
	// Sanger encoding is the default, but if solexaEncoding is true, we 
	// use the alternate version
	public static int fastqQualCodeToInt(char c, boolean solexaEncoding)
	{
		int v;
		if (solexaEncoding)
		{
			// Pre pipe 1.4
			// v = (int)Math.round(10.0 * Math.log10(1+Math.pow(10.0,(((int)c-64)/10.0))));
			 
			// Pipe >= 1.4
			v = (int)c-64;
		}
		else
		{
			v = (int)c - 33;
		}

		return v;
	}
		
	// Evaluate n!
    public static long factorial( int n )
    {
        if( n <= 1 )     // base case
            return 1;
        else
            return n * factorial( n - 1 );
    }


	public static boolean GffRecordsOverlap(GFFRecord a,	GFFRecord b) 
	{
		int aS = a.getStart();
		int bS = b.getStart();
		int aE = a.getEnd();
		int bE = b.getEnd();
		
		if ((aE<bS) || (bE<aS))
		{
			return false;
		}
		else
		{
			return true;
		}
	}
    
}
