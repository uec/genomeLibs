/**
 * 
 */
package edu.usc.epigenome.genomeLibs.ReadPos.Streamers;

import java.util.*;

import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;
import edu.usc.epigenome.genomeLibs.ReadPos.StreamHandlers.ReadPosStreamHandler;

/**
 * @author benb
 *
 */
public class ReadPosStreamer extends LinkedList<ReadPosStreamHandler>{

	/** Obj vars **/
	protected Iterator<ReadPos> rpIt;

	private static final long serialVersionUID = 7721L;
	
	/**
	 * @param inRpIt
	 */
	public ReadPosStreamer(Iterator<ReadPos> inRpIt)
	{
		rpIt = inRpIt;
	}
	
	
	/**
	 * 
	 */
	public void run()
	{
		// Call initialization functions of all handlers
		this.initHandlers();
		
		// Iterate over APs
		this.iterateRps();
		
		// Call finalizations functions of all handlers
		this.finishHandlers();
	}
	
	
	protected void initHandlers()
	{
		Iterator<ReadPosStreamHandler> handlerIt = this.iterator();
		while (handlerIt.hasNext())
		{
			handlerIt.next().init();
		}
	}
	
	
	protected void finishHandlers()
	{
		Iterator<ReadPosStreamHandler> handlerIt = this.iterator();
		while (handlerIt.hasNext())
		{
			handlerIt.next().finish();
		}
	}
	
	protected void iterateRps()
	{
		// Now go through one by one until we hit the end 
		while (rpIt.hasNext())
		{
			ReadPos curRp = rpIt.next();
			this.processRp(curRp);
		}
		
		// Finish up the current windows
	}
	
	protected boolean processRp(ReadPos currentAp)
	{
		boolean passes = true;

		Iterator<ReadPosStreamHandler> handlerIt = this.iterator();
		while (passes && handlerIt.hasNext())
		{
			passes &= handlerIt.next().streamElement(currentAp);
		}
		
		return passes;
	}

	
}
