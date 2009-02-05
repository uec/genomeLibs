/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers;

import java.util.*;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosNull;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.AlignmentPosStreamHandler;

/**
 * @author benb
 *
 */
public class AlignmentPosStreamer extends LinkedList<AlignmentPosStreamHandler> {

	/** Obj vars **/
	protected Iterator<AlignmentPos> apIt;

	protected int preWindSize;
	protected int postWindSize;
	private AlignmentPos[] apBuffer = null; // Temporary storage for the buffer
	private AlignmentPos[] preBuffer = null; // Temporary storage for the buffer
	private AlignmentPos[] postBuffer = null; // Temporary storage for the buffer
	private AlignmentPos[] preNull = null; // Temporary storage for the buffer
	private AlignmentPos[] postNull = null; // Temporary storage for the buffer
	
	
//	protected Queue queue = null;
	
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
		apBuffer = new AlignmentPos[preWindSize+postWindSize+1];
		preBuffer = new AlignmentPos[preWindSize];
		postBuffer = new AlignmentPos[postWindSize];
		preNull = new AlignmentPos[preWindSize];
		for (int i = 0; i < preWindSize; i++) { preNull[i] = new AlignmentPosNull(); }
		postNull = new AlignmentPos[postWindSize];
		for (int i = 0; i < postWindSize; i++) { postNull[i] = new AlignmentPosNull(); }
		
	}
	
	
	/**
	 * 
	 */
	
//	protected void initQueue()
//	{
//		// This could conceivable be either a LinkedList or Array based queue. Linked
//		// list probably has better throughput and is not a big deal for memory since
//		// the windows are generally small.
//		queue = new LinkedList(); 
//	}
	
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
	
	
	protected void finishHandlers()
	{
		Iterator<AlignmentPosStreamHandler> handlerIt = this.iterator();
		while (handlerIt.hasNext())
		{
			handlerIt.next().finish();
		}
	}
	
	protected void iterateAps()
	{
		// The general strategy is to start a new queue at the beginning 
		// of each chromosome.  When we start a new queue, we first fill it
		// up (including nulls to fill up the preWind slots), then we go until 
		// we get to the end of the chromosome or the last AP.  At the end, 
		// we add nulls to fill up the posWind buffer
		//
		// Note that the queue contains APs that are contiguous only in the
		// sense of being contiguous emissions by the iterator. They are not
		// necessarily in adjacent positions on the chromosome (they are
		// guaranteed to be in order however, as specified in the contract
		// of the AlignmentPosIterator). It is the job of the processQueue
		// function to make sure that arrays passed to the handlers contain
		// APs in adjacent chromosomal positions.
		
		String currentChr = AlignmentPos.NULL_CHROM;
		Queue<AlignmentPos> queue = null;
		
		while (apIt.hasNext())
		{
			AlignmentPos curAp = apIt.next();
			if (!curAp.getChr().equalsIgnoreCase(currentChr))
			{
				// System.err.println(curAp.getChr() + " != " + currentChr);
				
				// New chrom.  First, finish up the old chrom if necessary 
				if (queue != null)
				{
					finishChrom(queue);
				}
				
				// Now initialize and fill up the queue.
				// This could conceivable be either a LinkedList or Array based queue. Linked
				// list probably has better throughput and is not a big deal for memory since
				// the windows are generally small.
				queue = new LinkedList<AlignmentPos>();
				this.startChrom(queue, curAp);
			}
			else
			{
				// System.err.println(curAp.getChr() + " == " + currentChr);
				// Same chrom. Remove one from head, add one to tail, process.
				queue.poll();
				queue.add(curAp);
				processQueue(queue);
			}
			
			// Increment
			currentChr = curAp.getChr();
		}
		
		// Finish off our last chromosome
		finishChrom(queue);
	}
	
	
	/**
	 * Processes the first preWindSize elements at the beginning of a chromosome.
	 * May read from apIt and call processQueue
	 * 
	 * @param queue
	 * @param firstAp is optional.  If we have it, we use it and then take
	 *  the remaining (preWindSize-1) from apIt.  Otherwise, we take all 
	 *  preWindSize from apIt
	 */
	protected void startChrom(Queue<AlignmentPos> queue, AlignmentPos firstAp)
	{
//		System.err.println("startChrom(" + AlignmentPos.getRefTokens(queue) + ", " + ((firstAp==null) ? "null" : firstAp.getRefToken()) + ")");
//		System.err.println("startChrom(" + ((firstAp==null) ? "null" : firstAp.getChr()) + ")");

		// First push on the preWind null APs
		for (int i = 0; i < preWindSize; i++)
		{
			queue.add(new AlignmentPosNull());
		}
		
		// Now fill up the currentAp and postWind
		for (int i = 0 ; i < (postWindSize+1); i++)
		{
			if (firstAp != null)
			{
				queue.add(firstAp);
				firstAp = null;
			}
			else
			{
				AlignmentPos nextAp = (apIt.hasNext()) ? apIt.next() : (new AlignmentPosNull());
				queue.add(nextAp);
			}
		}
		
		// Process the queue
//		System.err.println("\tfinish startChrom(" + AlignmentPos.getRefTokens(queue) + ", " + ((firstAp==null) ? "null" : firstAp.getRefToken()) + ")");
		processQueue(queue);
	}
	
	/**
	 * Processes the remaining postWindSize elements at the end of a chromosome.
	 * Does not read from apIt but does call processQueue
	 * 
	 * @param queue
	 */
	protected void finishChrom(Queue<AlignmentPos> queue)
	{
//		System.err.println("finishChrom(" + AlignmentPos.getRefTokens(queue) + ")");
		
		for (int i = 0; i < postWindSize; i++)
		{
			queue.poll(); // Remove from head
			queue.add(new AlignmentPosNull());
			processQueue(queue);
		}
	}

	
	/**
	 * @param apQueue
	 * @return true if the element passes
	 */
	//TODO We could do a better job of managing buffers that are partially
	// contiguous
	protected boolean processQueue(Queue<AlignmentPos> apQueue)
	{
//		System.err.print(AlignmentPos.getRefTokens(apQueue) + "  ");

		
		// Copy the pre and post portions of the queue to separate arrays.
		apQueue.toArray(apBuffer);
		System.arraycopy(apBuffer, 0, preBuffer, 0, preWindSize);
		System.arraycopy(apBuffer, preWindSize+1, postBuffer, 0, postWindSize);
		AlignmentPos currentAp = apBuffer[preWindSize];
		
		AlignmentPos[] priorAps = preBuffer;
		AlignmentPos[] postAps = postBuffer;
		
		// Check if the preWindow and postWindow are actually in 
		// contiguous position in the alignment
		int currentApPos = currentAp.getPos();
		if ((preWindSize>0) && (priorAps[0].getPos() != (currentApPos-preWindSize)))
		{
//			System.err.println("PreWind non-contiguous: priorAps= " +  priorAps[0].getPos() +
//					", currentPos=" + currentApPos + ", preWindSize=" + preWindSize);
			priorAps = preNull;
		}
		if ((postWindSize>0) && (postAps[postWindSize-1].getPos() != (currentApPos+postWindSize)))
		{
//			System.err.println("PostWind non-contiguous " + currentAp.toString());
			postAps = postNull;
		}
		
		return processAp(priorAps, currentAp, postAps);
	}
	

	
	
	/**
	 * @param priorAps
	 * @param currentAp
	 * @param postAps
	 * @return
	 * 
	 * priorAps and nextAps are guaranteed to contain APs which are adjacent to each other 
	 * on the chromosome.
	 */
	protected boolean processAp(AlignmentPos[] priorAps, 
			AlignmentPos currentAp, AlignmentPos[] postAps)
	{
		// Now pass it through handler
		boolean passes = true;
		Iterator<AlignmentPosStreamHandler> handlerIt = this.iterator();
		while (passes && handlerIt.hasNext())
		{
			AlignmentPosStreamHandler handler = handlerIt.next();
////			System.err.println("Streaming to " + handler + ":\t" + 
//			System.err.println(AlignmentPos.getRefTokens(priorAps) + 
//					"," + currentAp.getRefToken() + "," + AlignmentPos.getRefTokens(postAps));
			passes &= handler.streamElement(priorAps, currentAp, postAps);
		}
		return passes;
	}

	
}
