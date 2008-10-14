package edu.usc.epigenome.genomeLibs;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;

public abstract class AlignmentPosIterator implements Iterator<AlignmentPos> {

	protected BufferedReader f_open_stream = null;
	protected String f_open_file = null;
	protected AlignmentPosOptions f_options = null;
	
	public AlignmentPosIterator(String fn, AlignmentPosOptions apos)
	throws IOException
	{
		f_open_stream = new BufferedReader(new FileReader(fn));
		f_open_file = fn;
		f_options = apos;
	}
	
	/***********
	 * Interface functions
	 */
	
	public boolean hasNext()
	{
		boolean out = false;
		
		try
		{
			out = (f_open_stream != null) && (f_open_stream.ready());
		}
		catch (Exception e)
		{
			System.err.println("Could not read file " + f_open_file + "\n" + e.toString());
		}
		
		return out;
	}
	
	public AlignmentPos next()
	{
		AlignmentPos ap = null; 
		
		try
		{
			ap = nextAlignment();
		}
		catch (Exception e)
		{
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

	abstract protected AlignmentPos nextAlignment() throws IOException;
	
}
