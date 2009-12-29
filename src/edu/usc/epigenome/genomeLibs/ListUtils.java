package edu.usc.epigenome.genomeLibs;

import java.util.*;
import java.io.*;
import java.nio.*;

public class ListUtils {

	protected static String delim = ",";

	//TODO Not thread safe
	protected static final int STRINGBUFLEN = 100000000;
	protected static StringBuilder STRINGBUF = new StringBuilder(STRINGBUFLEN);

	public static void setDelim(String d)
	{
		delim = d;
	}

	
	//TODO Not threadsafe because of static buffer STRINGBUF
	public static String[] readLineSplitByChar(Reader r, char delim, int countEstimate)
	throws Exception
	{
		Vector<String> v = new Vector<String>(countEstimate);

		STRINGBUF.delete(0, STRINGBUFLEN);
		int c;

		int onC = 0;
		
		charLoop:
		while ((c = r.read()) >= 0)
		{
			if (c == delim)
			{
				v.add(new String(STRINGBUF));
				STRINGBUF.delete(0, STRINGBUFLEN);
				onC=0;
			}
			else if ((c == '\n') || (c == '\r')) 
			{
				break charLoop;
			}
			else
			{
				STRINGBUF.append((char)c);
			}
			onC++;
			if (onC>=STRINGBUFLEN)
				throw new Exception("Can not process lines with a field more than " + STRINGBUFLEN + " characters");
		}
		
		// Add the last one (if there is one)
		if (STRINGBUF.length()>0) v.add(new String(STRINGBUF));
		
		return (String[])v.toArray( new String[v.size()] );
	}
	
	
	
	// Is Java split inefficient??
	public static String[] splitByChar(String in, char delim)
	{
		char[] ar = in.toCharArray();
		int len = ar.length;
		
		Vector<String> vec = new Vector<String>(1000);
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
		StringBuffer sb = new StringBuffer(l.length*5);
		//System.err.println("excelLine length = " + l.length);
		
		for (int i = 0; i < l.length; i++)
		{
			if (i > 0) sb.append(delim);
			//System.err.println("\tOutputting " + l[i]);
			//sb.append(String.format("%g", l[i]));
			sb.append(l[i]);
		}
		//System.err.println(sb.toString());
		return sb.toString();
	}
	
	public static String excelLine(int[] l)
	{
		STRINGBUF.setLength(0);
		for (int i = 0; i < l.length; i++)
		{
			if (i > 0) STRINGBUF.append(delim);
			STRINGBUF.append(l[i]);
		}
		return STRINGBUF.toString();
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
		
	public static String excelLine(List<Object> l)
	{
		ListIterator<Object> li = l.listIterator();
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
	
	public static Vector<String[]> parseExcelFile(String fn)
	throws Exception
	{
		BufferedReader in = new BufferedReader(new FileReader(fn));
		
		Vector<String[]> out = new Vector<String[]>();
		String line;
		while ((line=in.readLine()) != null)
		{
			out.add(parseExcelLine(line));
		}
		
		return out;
	}
	
	public static Vector<String> parseSingleColumnFile(String fn, boolean header)
	throws Exception
	{
		return parseExcelFileColumn(fn, 1, header);
	}
	
	
	public static Vector<String> parseExcelFileColumn(String fn, int col, boolean header)
	throws Exception
	{
		Vector<String[]> lines = parseExcelFile(fn);
		
		Vector<String> out = new Vector<String>();
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

	public static String tabbedLine(List<Object> l)
	{
		ListIterator<Object> li = l.listIterator();
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
	
	
	public static double[] intArrToDouble(int[] in)
	{
		double[] out = new double[in.length];
		
		for (int i = 0; i < in.length; i++)
		{
			out[i] = (double)in[i];
		}
		
		return out;
	}
	                  
	                  

}
