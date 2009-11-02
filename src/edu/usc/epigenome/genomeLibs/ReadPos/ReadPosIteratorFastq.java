package edu.usc.epigenome.genomeLibs.ReadPos;

import java.io.IOException;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

import edu.usc.epigenome.genomeLibs.MiscUtils;

/**
 * @author benb
 *
 * These are fake alignment positions, actually just individual read positions
 *
 */


public class ReadPosIteratorFastq extends ReadPosIterator {

	//private int totalBasesRead = 0;
	static final protected int IN_NO_SEC = 0;
	static final protected int IN_NUC_SEC = 1;
	static final protected int IN_QUAL_SEC = 2;
	static final protected int READ_BUFFERED = 3;
	static final protected int END_OF_FILE = 4;

	// state vars
	protected String currentNucs = null;
	protected int currentNucsCycle = 0;
	protected String currentQuals = null;
	
	
	public ReadPosIteratorFastq(String fn, ReadPosOptions apos) 
	throws IOException {
		super(fn, apos);
	}

	protected void initStoredNucs()
	{
		currentNucs = null;
		currentNucsCycle = 0;
		currentQuals = null;
	}
	
	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.ReadPosIterator#hasNext()
	 */
	@Override
	public boolean hasNext() 
	{
		
		boolean out = false;
		try
		{
			ReadPos next = this.nextReadPos(true);
			out = (next != null);
		}
		catch (Exception e)
		{
			System.err.println("Could not get next position");
			e.printStackTrace();
			System.exit(0);
		}
		
		//System.err.println("Has next? " + out);
		return out;
	}

	@Override
	protected ReadPos nextReadPos()
	throws IOException, IllegalSymbolException, Exception
	{
		//System.err.println("Calling nextReadPos()");
		return nextReadPos(false);
	}

	protected ReadPos nextReadPos(boolean rollback)
	throws IOException, IllegalSymbolException, Exception
	{
		int state = IN_NO_SEC;
		ReadPos outRp = null;

		while ((outRp == null) && (state != END_OF_FILE))
		{
			// We'll check if we have already started one sequence read.  If not,
			// we will go fetch the next one.
			outRp = fetchNextReadPosStored(rollback);
			if (outRp==null)
			{
				// If we don't have a stored one, go until we fetch the next one
				this.initStoredNucs();
				state = IN_NO_SEC;
				while ((state != READ_BUFFERED) && (state != END_OF_FILE))
				{
					String line = this.openStream.readLine();
					//System.err.println("Read line: " + line);
					if (line == null)
					{
						state = END_OF_FILE;
					}
					else if (line.length() == 0)
					{
						// do nothing
					}
					else
					{
						char firstChar = line.charAt(0);

						switch (firstChar)
						{
						case '@':
							state = IN_NUC_SEC;
							break;
						case '#':
							// A comment, do nothing
							break;
						case '+':
							state = IN_QUAL_SEC;
							break;
						default:
							if (state == IN_NUC_SEC)
							{
								currentNucs = line;
								// After the nuc sec, we go back to empty state
								state = IN_NO_SEC;
							}
							else if (state == IN_QUAL_SEC)
							{
								currentQuals = line;
								// After the qual sec, we are done
								state = READ_BUFFERED;
							}
							else
							{
								throw new Exception("Error parsing fastq file: Found an illegal line " + 
										"which does not start with a @ or + character:\n" + line);
							}
						break;
						}
					}
					
				}
			}
		}
		
//		System.err.println("outRp = " + ((outRp==null) ? "null" : outRp.commaSeparatedLine()) + "\tstate="+ state + "\trollback=" + rollback);
		return outRp;
	}
	

	/**
	 * If there is a currentNucs/currentQuals/currentNucsCycle state, we take the
	 * next nuc/qual pair based on cycle number, then increment currentNucsCycle. 
	 * 
	 * @rollback set to true if you don't want to advance to next one
	 * 
	 * @return ReadPos, or null if there is not another nucleotide stored
	 */
	protected ReadPos fetchNextReadPosStored(boolean rollback)
	throws IllegalSymbolException
	{
		ReadPos outRp = null;
		if (currentNucs != null)
		{
			// We have to loop because we might have read positions with 
			// bad quality scores.
			boolean done = false;
			while (!done)
			{

				// Pop one off and serve it (unless we hit the end, in which case we
				// just initialize.
				if (currentNucsCycle >= currentNucs.length())
				{
					this.initStoredNucs();
					done = true;
				}
				else
				{
					Symbol sym = DNATools.forSymbol(currentNucs.charAt(currentNucsCycle));
					int qual = MiscUtils.fastqQualCodeToInt(currentQuals.charAt(currentNucsCycle),this.rpOptions.positionQualsSolexaEncoding);
					//System.err.println("\t" + currentNucs.charAt(currentNucsCycle) + "\t" + qual + "  (" + currentQuals.charAt(currentNucsCycle) + ") >= " + rpOptions.minQualityScore);
					if (qual >= rpOptions.minQualityScore)
					{
						outRp = new ReadPos(sym,StrandedFeature.UNKNOWN);
						if (rpOptions.trackPositions || rpOptions.trackQuals)
						{
							// add 1 to currentNucsCycle to go from 0-based to 1-based
							int rpCycle = (rpOptions.trackPositions) ? (currentNucsCycle+1) : ReadPos.UNKNOWN_CYCLE;
							int rpQual = (rpOptions.trackQuals) ? qual : ReadPos.UNKNOWN_QUAL;
							outRp = new ReadPosRich(outRp, rpCycle, rpQual);
						}
						if (rollback) currentNucsCycle--;
					}
					done = (outRp != null);

					// And increment to the next nuc
					currentNucsCycle++;
				}
			}
		}

		return outRp;
	}
	
//	protected boolean haveStoredReadPos()
//	{
//		boolean out = false;
//		
//		try
//		{
//			ReadPos rp = fetchNextReadPosStored(true);
//			out = (rp != null);
//		}
//		catch (IllegalSymbolException e)
//		{
//			System.err.println("IllegalSymbolException trying to read next line\n");
//			e.printStackTrace();
//			System.exit(0);
//		}
//
//		return out;
//	}
	
	
}
