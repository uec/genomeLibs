package edu.usc.epigenome.genomeLibs.ReadPos;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;

public abstract class ReadPosIterator implements Iterator<ReadPos> {

	protected BufferedReader openStream = null;
	protected String openFile = null;
	protected ReadPosOptions rpOptions = null;
	protected int numRead = 0;
	
	
	public ReadPosIterator(String fn, ReadPosOptions rpos)
	throws IOException
	{
		openStream = new BufferedReader(new FileReader(fn));
		openFile = fn;
		rpOptions = rpos;
	}
	
	/***********
	 * Interface functions
	 */
	
	abstract public boolean hasNext();

	
	public ReadPos next()
	{
		ReadPos rp = null; 
		
		try
		{
			rp = nextReadPos();
			this.numRead++;
			if ((this.numRead % 1000000) == 0) System.err.println(this.numRead + " APs processed");
		}
		catch (Exception e)
		{
			System.err.println("Could not get next ReadPos:\n");
			e.printStackTrace();
			System.exit(0);
		}
		
		return rp;
	}
	
	public void remove()
	throws UnsupportedOperationException 
	{
		throw new UnsupportedOperationException();
	}
	
	/**********
	 * Unimplemented functions
	 */

	abstract protected ReadPos nextReadPos() throws Exception;
	
}
