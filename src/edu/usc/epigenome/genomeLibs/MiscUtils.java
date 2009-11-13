package edu.usc.epigenome.genomeLibs;


public class MiscUtils {


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
		

}
