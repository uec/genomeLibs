package edu.usc.epigenome.genomeLibs;

import java.util.*;
import java.io.*;

public class ListUtils {

	protected static String delim = ",";

	public static void setDelim(String d)
	{
		delim = d;
	}

	
	// Is Java split inefficient??
	public static String[] splitByChar(String in, char delim)
	{
		char[] ar = in.toCharArray();
		int len = ar.length;
		
		Vector vec = new Vector(1000);
		char[] buf = new char[len];
		int buf_pos = 0;
		for (int i = 0; i < len; i++)
		{
			char c = ar[i];
			
			if ( c ==  delim )
			{
				vec.add(new String(buf,0,buf_pos));
				buf_pos = 0;
			}
			else
			{
				buf[buf_pos++] = c;
			}
		}
		
		// Add the last one
		if (buf_pos>0) vec.add(new String(buf,0,buf_pos));

		return (String[])vec.toArray( new String[vec.size()] );
	}
	
	public static String excelLine(String[] l)
	{

		String out = "";
		for (int i = 0; i < l.length; i++)
		{
			if (i > 0) out += delim;
			out += (l[i] == null) ? "" : l[i];
		}
		return out;
	}
	
	public static String excelLine(double[] l)
	{

		String out = "";
		for (int i = 0; i < l.length; i++)
		{
			if (i > 0) out += delim;
			out += l[i];
		}
		return out;
	}
	
	public static String excelLine(int[] l)
	{

		String out = "";
		for (int i = 0; i < l.length; i++)
		{
			if (i > 0) out += delim;
			out += l[i];
		}
		return out;
	}
	
	public static String excelLine(boolean[] l)
	{

		String out = "";
		for (int i = 0; i < l.length; i++)
		{
			// if (i > 0) out += delim;
			out +=  (l[i]) ? 1 : 0;
		}
		return out;
	}
	
	public static String excelLine(List l)
	{
		ListIterator li = l.listIterator();
		String out = "";
		while (li.hasNext())
		{
			out += li.next();
			if (li.hasNext())
			{
				out += delim;
			}
		}
		return out;
	}
	
	public static Vector parseExcelFile(String fn)
	throws Exception
	{
		BufferedReader in = new BufferedReader(new FileReader(fn));
		
		Vector out = new Vector();
		String line;
		while ((line=in.readLine()) != null)
		{
			out.add(parseExcelLine(line));
		}
		
		return out;
	}
	
	public static Vector parseSingleColumnFile(String fn, boolean header)
	throws Exception
	{
		return parseExcelFileColumn(fn, 1, header);
	}
	
	
	public static Vector parseExcelFileColumn(String fn, int col, boolean header)
	throws Exception
	{
		Vector lines = parseExcelFile(fn);
		
		Vector out = new Vector();
		int i = (header) ? 1 : 0;
		while (i < lines.size())
		{
			String[] line = (String[])lines.get(i);
			out.add(line[col-1]);
			i++;
		}
		return out;
	}
	
	public static String[] parseExcelLine(String s)
	{
		String delim_full = delim;
		if (delim.equals(","))
		{
			delim_full = "[ \t]*" + delim + "[ \t]*";
		}
		String[] result = s.split(delim_full);
		return result;
	}

	public static String tabbedLine(List l)
	{
		ListIterator li = l.listIterator();
		String out = "";
		while (li.hasNext())
		{
			out += li.next();
			if (li.hasNext())
			{
				out += "\t";
			}
		}
		return out;
	}

}
