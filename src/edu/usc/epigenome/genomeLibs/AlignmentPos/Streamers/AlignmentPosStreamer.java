/**
 * 
 */
package edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers;

import java.util.*;

import org.biojava.bio.symbol.Symbol;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosNull;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.AlignmentPosStreamHandler;
import edu.usc.epigenome.genomeLibs.Counters.NmerCounter;
import edu.usc.epigenome.genomeLibs.Counters.StringCounter;
import edu.usc.epigenome.genomeLibs.Counters.SymbolCounter;

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
	
	private NmerCounter preNmerCounter = null; // Keep track of symbol counts
	private NmerCounter postNmerCounter = null; // Keep track of symbol counts
	private boolean nmerCounterStale = true;
	private Symbol[] singleSymBuf = new Symbol[1];
	private Symbol[] doubleSymBuf = new Symbol[2];
	
	protected AlignmentPos sApNull = new AlignmentPosNull(); // so that we don't have to create a new one each time.
	
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
		for (int i = 0; i < preWindSize; i++) { preNull[i] = sApNull; }
		postNull = new AlignmentPos[postWindSize];
		for (int i = 0; i < postWindSize; i++) { postNull[i] = sApNull; }
		
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
				AlignmentPos dropped = queue.poll();
				queue.add(curAp);
				processQueue(queue, dropped);
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
				AlignmentPos nextAp = (apIt.hasNext()) ? apIt.next() : sApNull;
				queue.add(nextAp);
			}
		}
		
		// Process the queue
		// System.err.println("\tfinish startChrom(" + AlignmentPos.getRefTokens(queue) + ", " + ((firstAp==null) ? "null" : firstAp.getRefToken()) + ")");
		symbolCountersMakeStale();
		processQueue(queue, null);
	}
	
	/**
	 * Processes the remaining postWindSize elements at the end of a chromosome.
	 * Does not read from apIt but does call processQueue
	 * 
	 * @param queue
	 */
	protected void finishChrom(Queue<AlignmentPos> queue)
	{
		// System.err.println("finishChrom(" + AlignmentPos.getRefTokens(queue) + ")");
		
		for (int i = 0; i < postWindSize; i++)
		{
			AlignmentPos dropped = queue.poll(); // Remove from head
			queue.add(new AlignmentPosNull());
			processQueue(queue, dropped);
		}
	}

	
	/**
	 * @param apQueue
	 * @return true if the element passes
	 */
	//TODO We could do a better job of managing buffers that are partially
	// contiguous
	protected boolean processQueue(Queue<AlignmentPos> apQueue, AlignmentPos droppedAp)
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
		int windMinPos = currentApPos-preWindSize;
		int windMaxPos = currentApPos+postWindSize;
		
		// Invalidate APs that might be in the pre window buffer but not actually in the
		// window (due to non-contiguity)
		if ((preWindSize>0) && (priorAps[0].getPos() != windMinPos))
		{
//			System.err.println("PreWind non-contiguous: priorAps= " +  priorAps[0].getPos() +
//					", currentPos=" + currentApPos + ", preWindSize=" + preWindSize);
			
			boolean done = false;
			for (int ind=0; !done && (ind<priorAps.length); ind++)
			{
				AlignmentPos ap = priorAps[ind];
				if (ap.getPos() == -1) // AlignmentPosNull, skip it
				{
					
				}
				else if (!ap.getChr().equals(currentAp.getChr()) || (ap.getPos() < windMinPos))
				{
					priorAps[ind] = sApNull;
				}
				else	
				{
					done = true;
				}
			}

			this.symbolCountersMakeStale();
		}

		// Invalidate APs that might be in the post window buffer but not actually in the
		// window (due to non-contiguity)
		if ((postWindSize>0) && (postAps[postWindSize-1].getPos() != windMaxPos))
		{
			// System.err.println("PostWind non-contiguous " + currentAp.toString());

			boolean done = false;
			for (int ind=(postAps.length)-1; !done && (ind>=0); ind--)
			{
				AlignmentPos ap = postAps[ind];
				//System.err.println("\tTesting #" + ind + " , " + ap.getPos());
				if (ap.getPos() == -1) // AlignmentPosNull, skip it
				{
				}
				else if (!ap.getChr().equals(currentAp.getChr()) || (ap.getPos() > windMaxPos))
				{
					postAps[ind] = sApNull;
				}
				else	
				{
					done = true;
				}
			}		
			
			this.symbolCountersMakeStale();
		}
		
		// Pass to handlers
		AlignmentPosStreamerPosition streamPos = new AlignmentPosStreamerPosition(preWindSize,postWindSize);
		streamPos.priorAps = priorAps;
		streamPos.currentAp = currentAp;
		streamPos.nextAps = postAps;
		
		// Subclasses or handlers might want to override this.  But by default it's the AP itself.
		streamPos.currentScoredPosition = currentAp;
		
		// Increment and pass on the nmer counters.
		// This is VERY slow, due to NmerCounter.keyFromSyms
		if (false)
		{
			this.symbolCountersIncrement(droppedAp, priorAps, currentAp, postAps);
			streamPos.preNmerCounts = (this.symbolCountersAreStale()) ? null : this.preNmerCounter;
			streamPos.nextNmerCounts = (this.symbolCountersAreStale()) ? null : this.postNmerCounter;
		}
		
		return processAp(streamPos);
	}
	

	
	
	/**
	 * @param streamPos
	 * @return
	 * 
	 * priorAps and nextAps are guaranteed to contain APs which are adjacent to each other 
	 * on the chromosome.
	 */
	protected boolean processAp(AlignmentPosStreamerPosition streamPos)
	{
		// Now pass it through handler
		boolean passes = true;
		Iterator<AlignmentPosStreamHandler> handlerIt = this.iterator();
		while (passes && handlerIt.hasNext())
		{
			AlignmentPosStreamHandler handler = handlerIt.next();
//			// System.err.print("Streaming to " + handler + ":\t"); 
//			System.err.println(AlignmentPos.getRefTokens(streamPos.priorAps) + 
//					"," + streamPos.currentAp.getRefToken() + "," + AlignmentPos.getRefTokens(streamPos.nextAps));
			passes &= handler.streamElement(streamPos);
		}
		return passes;
	}

	
	
	
	// -------------------------------  Symbol count handling ------------------
	
	protected void symbolCountersMakeStale()
	{
		nmerCounterStale = true;
	}
	
	protected void symbolCountersMakeFresh()
	{
		nmerCounterStale = false;
	}

	protected boolean symbolCountersAreStale()
	{
		return nmerCounterStale;
	}
	
	protected void symbolCountersInit(AlignmentPos[] preAps, AlignmentPos currentAp, AlignmentPos nextAps[])
	{
		if ((preAps != null) && (preWindSize>0))
		{
			preNmerCounter = new NmerCounter();
			preNmerCounter.incrementAllSubstrings(AlignmentPos.getSymbols(preAps),2);
		}

		if ((nextAps != null) && (postWindSize > 0))
		{
			postNmerCounter = new NmerCounter();
			postNmerCounter.incrementAllSubstrings(AlignmentPos.getSymbols(nextAps),2);
		}
	}

	protected void symbolCountersIncrement(AlignmentPos droppedAp, AlignmentPos[] preAps, AlignmentPos currentAp, AlignmentPos nextAps[])
	{
		if (symbolCountersAreStale())
		{
			symbolCountersInit(preAps, currentAp, nextAps);
			this.symbolCountersMakeFresh();
		}
		else
		{
			if ((preAps!=null) && (preWindSize>0))
			{
				// Increment pre
				singleSymBuf[0] = preAps[preAps.length-1].getRef();
				preNmerCounter.increment(singleSymBuf);
				if (preWindSize>1)
				{
					doubleSymBuf[0] = preAps[preAps.length-2].getRef();
					doubleSymBuf[1] = preAps[preAps.length-1].getRef();
					preNmerCounter.increment(doubleSymBuf);
				}

				// Decrement pre
				singleSymBuf[0] = droppedAp.getRef();
				preNmerCounter.decrement(singleSymBuf);
				if (preWindSize>1)
				{
					doubleSymBuf[0] = droppedAp.getRef();
					doubleSymBuf[1] = preAps[0].getRef();
					preNmerCounter.decrement(doubleSymBuf);
				}
			}

			if ((nextAps != null) && (postWindSize>0))
			{
				// Increment post
				singleSymBuf[0] = nextAps[nextAps.length-1].getRef();
				postNmerCounter.increment(singleSymBuf);
				if (postWindSize>1)
				{
					doubleSymBuf[0] = nextAps[nextAps.length-2].getRef();
					doubleSymBuf[1] = nextAps[nextAps.length-1].getRef();
					postNmerCounter.increment(doubleSymBuf);
				}

				// Decrement post
				singleSymBuf[0] = currentAp.getRef();
				postNmerCounter.decrement(singleSymBuf);
				if (postWindSize>1)
				{
					doubleSymBuf[0] = currentAp.getRef();
					doubleSymBuf[1] = nextAps[0].getRef();
					postNmerCounter.decrement(doubleSymBuf);
				}
			}
		}
	}

	
}
