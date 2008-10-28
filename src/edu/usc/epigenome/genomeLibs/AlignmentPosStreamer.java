/**
 * 
 */
package edu.usc.epigenome.genomeLibs;

import java.util.*;

/**
 * @author benb
 *
 */
public class AlignmentPosStreamer extends LinkedList<AlignmentPosStreamHandler>{

	/** Obj vars **/
	protected Iterator<AlignmentPos> apIt;

	protected int preWindSize;
	protected int postWindSize;
	protected LinkedList<AlignmentPos> preWindow = new LinkedList<AlignmentPos>(); // Linked list for forwared and reverse navigation
	protected LinkedList<AlignmentPos> postWindow = new LinkedList<AlignmentPos>(); // Linked list for forwared and reverse navigation
	
	private static final long serialVersionUID = 7720L;
	
	/**
	 * @param inApIt
	 * @param inPreWindSize
	 * @param inPostWindSize
	 */
	public AlignmentPosStreamer(Iterator<AlignmentPos> inApIt, int inPreWindSize, int inPostWindSize)
	{
		apIt = inApIt;
		preWindSize = inPreWindSize;
		postWindSize = inPostWindSize;
	}
	
	
	/**
	 * 
	 */
	public void run()
	{
		// Call initialization functions of all handlers
		this.initHandlers();
		
		// Iterate over APs
		this.iterateAps();
		
		// Call finalizations functions of all handlers
		this.finishHandlers();
	}
	
	
	protected void initHandlers()
	{
		Iterator<AlignmentPosStreamHandler> handlerIt = this.iterator();
		while (handlerIt.hasNext())
		{
			handlerIt.next().init();
		}
	}
	
	protected void iterateAps()
	{
		// The "next" AP is actually the end of the post window.  We must do window
		// adjustment at the beginning of each new chromosome
		String currentChr = "noChrom";
		
		// Now go through one by one until we hit the end 
		while (apIt.hasNext())
		{
			AlignmentPos curAp = apIt.next();
		
			// Are we on a new chrom?  If so, 
			
		}
		
		// Finish up the current windows
	}
	
	protected boolean processAp(LinkedList<AlignmentPos> priorAps, 
			AlignmentPos currentAp, LinkedList<AlignmentPos> nextAps)
	{
		boolean passes = true;

		Iterator<AlignmentPosStreamHandler> handlerIt = this.iterator();
		while (passes && handlerIt.hasNext())
		{
			passes &= handlerIt.next().streamElement(priorAps, currentAp, nextAps);
		}
		
		return passes;
	}
	
	protected void finishHandlers()
	{
		Iterator<AlignmentPosStreamHandler> handlerIt = this.iterator();
		while (handlerIt.hasNext())
		{
			handlerIt.next().finish();
		}
	}
	
}
