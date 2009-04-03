package edu.usc.epigenome.genomeLibs.AlignmentPos;

import java.io.IOException;

import org.biojava.bio.symbol.IllegalSymbolException;


import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.io.SeqIOTools;

public class AlignmentPosIteratorFasta extends AlignmentPosIterator {

	//private int totalBasesRead = 0;
	
	private SequenceIterator seqIt = null;
	
	private String curSeqResidues = null;
	private int curIdx = 0;
	private String curSeq = null;
	
	public AlignmentPosIteratorFasta(String fn) 
	throws IOException {
		super(fn, new AlignmentPosOptions());
		
		seqIt = SeqIOTools.readFastaDNA(this.openStream);
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosIterator#hasNext()
	 */
	@Override
	public boolean hasNext() {
		
		return (((curSeqResidues != null) && (curSeqResidues.length() > curIdx)) ||
				seqIt.hasNext());
	}

	@Override
	protected AlignmentPos nextAlignment()
	throws IOException, IllegalSymbolException, Exception
	{
		return nextAlignment(false);
	}

	protected AlignmentPos nextAlignment(boolean rollback)
	throws IOException, IllegalSymbolException, Exception
	{
		AlignmentPos ap = null;
		
		// If we're out of residues on the current sequence, get a new one.
		if ((curSeqResidues == null) || (curSeqResidues.length() <= curIdx))
		{
			Sequence seq = seqIt.nextSequence();
			if (seq==null) return null;
			curSeqResidues = seq.seqString();
			curSeq = seq.getName();
			curIdx = 0;
			// System.err.println("New seq " + curSeq);
		}
		
		if (curSeqResidues.length() <= curIdx) return null;

		char ref = curSeqResidues.charAt(curIdx);
		ap = new AlignmentPosRefOnly(ref, curSeq, curIdx+1);
		curIdx++;

		
		if (rollback) curIdx--;
		return ap;
	}
	
	
	

}
