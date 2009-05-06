package edu.usc.epigenome.genomeLibs.AlignmentPos;

import java.io.*;
import java.util.zip.*;

import java.util.GregorianCalendar;
import java.util.Iterator;

import org.jfree.util.StringUtils;


public abstract class AlignmentPosIterator implements Iterator<AlignmentPos> {

	protected BufferedReader openStream = null;
	protected String openFile = null;
	protected AlignmentPosOptions apOptions = null;
	protected int numRead = 0;
	
	
	public AlignmentPosIterator(String fn, AlignmentPosOptions apos)
	throws IOException
	{
		if (fn.endsWith(".gz"))
		{
			System.err.println("UNzipping " + fn);
			InputStream is = new GZIPInputStream(new FileInputStream(fn));
			InputStreamReader isr = new InputStreamReader(is);
			openStream = new BufferedReader(isr);
		}
		else
		{
			openStream = new BufferedReader(new FileReader(fn));
		}	
		
		openFile = fn;
		apOptions = apos;
	}
	
	/***********
	 * Interface functions
	 */
	
	public boolean hasNext()
	{
		boolean out = false;
		
		try
		{
			out = (openStream != null) && (openStream.ready());
		}
		catch (Exception e)
		{
			System.err.println("Could not read file " + openFile + "\n" + e.toString());
			System.exit(0);
		}
		
		return out;
	}
	
	public AlignmentPos next()
	{
		AlignmentPos ap = null; 
		
		try
		{
			ap = nextAlignment();
			this.numRead++;
			if ((this.numRead % 1000000) == 0)
//			if ((this.numRead % 1000) == 0)
			{
				System.err.println((new GregorianCalendar()).getTime() + "\t" + this.numRead + " APs processed");
			}
		}
		catch (Exception e)
		{
			System.err.println("Could not get next AlignmentPos:\n");
			e.printStackTrace();
			System.exit(0);
		}
		
		return ap;
	}
	
	public void remove()
	throws UnsupportedOperationException 
	{
		throw new UnsupportedOperationException();
	}
	
	/**********
	 * Unimplemented functions
	 */

	abstract protected AlignmentPos nextAlignment() throws Exception;
	
}
