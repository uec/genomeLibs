package edu.usc.epigenome.genomeLibs.AlignmentPos;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;


public abstract class AlignmentPosIterator implements Iterator<AlignmentPos> {

	protected BufferedReader openStream = null;
	protected String openFile = null;
	protected AlignmentPosOptions apOptions = null;
	protected int numRead = 0;
	
	
	public AlignmentPosIterator(String fn, AlignmentPosOptions apos)
	throws IOException
	{
		openStream = new BufferedReader(new FileReader(fn));
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
			if ((this.numRead % 1000000) == 0) System.err.println(this.numRead + " APs processed");
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
