package edu.usc.epigenome.genomeLibs.ChromScores;

import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.GoldAssembly;
import edu.usc.epigenome.genomeLibs.MiscUtils;

public class ChromScoresMotifPositions extends ChromScoresArrayInt {

	public ChromScoresMotifPositions(String genome) 
	throws Exception
	{
		super(genome);
		System.err.println("Initing GCConent with genome=" + genome);
	}

	
	public void populate(String target_chr, String motif, StrandedFeature.Strand desiredStrand)
	throws Exception
	{
		populate(target_chr, motif, 0, 0, desiredStrand);
	}

	
	
	// If end_coord == 0 , we do the whole chromosome 
	public void populate(String target_chr, String motif, int start_coord, int end_coord, StrandedFeature.Strand desiredStrand)
	throws Exception
	{
		// Get the chromosome
		char[] full_residues_arr = null;
		{
			Sequence seq = GoldAssembly.chromSeq(f_genome, target_chr);
			String full_residues = seq.seqString().toUpperCase();
			if (end_coord > 0)
			{
				full_residues = full_residues.substring(start_coord, end_coord);
			}	
			full_residues_arr = full_residues.toCharArray();
		}

		populate(target_chr, motif, full_residues_arr, start_coord, desiredStrand);
	}
	
	public void populate(String target_chr, String motif, char[] full_residues_arr, int start_coord, StrandedFeature.Strand desiredStrand)
	throws Exception
	{
		motif = motif.toUpperCase();
		char[] motifArr = motif.toCharArray();
		int motifLen = motif.length();
		int lastPos = full_residues_arr.length - motifLen;
		Integer oneAsInt = new Integer(1);
		
		boolean doFw = (desiredStrand == StrandedFeature.POSITIVE) || (desiredStrand == StrandedFeature.UNKNOWN);
		boolean doRev = (desiredStrand == StrandedFeature.NEGATIVE) || (desiredStrand == StrandedFeature.UNKNOWN);
		
		char[] motifRevcompArr = MiscUtils.revCompNucStr(motif).toCharArray();
		
		System.err.printf("LastPos=%d\n",lastPos);
		for (int startPos = 0; startPos<=lastPos; startPos++)
		{
			if (doFw)
			{
				boolean match = true;
				for (int motifPos = 0; match && (motifPos<motifLen); motifPos++)
				{
					int refPos = motifPos + startPos;
					match = nucMatches(motifArr[motifPos],full_residues_arr[refPos]);
					//match = (full_residues_arr[refPos] == motifArr[motifPos]);
				}
				if (match)
				{
					//System.err.printf("Found a (+) strand match to %s at pos %d\n", motif, startPos+start_coord);
					this.addScore(target_chr, startPos+start_coord, oneAsInt);
				}
			}
			
			if (doRev)
			{
				boolean match = true;
				for (int motifPos = 0; match && (motifPos<motifLen); motifPos++)
				{
					int refPos = motifPos + startPos;
					match = nucMatches(motifRevcompArr[motifPos],full_residues_arr[refPos]);
					//match = (full_residues_arr[refPos] == motifRevcompArr[motifPos]);
				}
				if (match)
				{
					//System.err.printf("Found a (-) strand match to %s at pos %d\n", motif, start_coord+startPos+motifLen-1);
					// This is a real problem.  Some of the MNase libraries seem to match start+len-1 (like?) start+len+1 (like ??), and others start+len (like Schones) 
					this.addScore(target_chr, start_coord+startPos+motifLen, oneAsInt);
				}
			}
		}
	}
	
	public static boolean nucMatches(char matchNuc, char actualNuc)
	throws Exception
	{
		boolean matches = false;
		
		switch(matchNuc)
		{
		case 'A':
		case 'C':
		case 'T':
		case 'G':
			matches = (actualNuc == matchNuc);
			break;
		case 'Y':
			switch(actualNuc)
			{
			case 'T':
			case 'C':
				matches = true;
				break;
			}
			break;
		case 'R':
			switch(actualNuc)
			{
			case 'A':
			case 'G':
				matches = true;
				break;
			}
			break;
		case 'W':
			switch(actualNuc)
			{
			case 'A':
			case 'T':
				matches = true;
				break;
			}
			break;
		case 'S':
			switch(actualNuc)
			{
			case 'C':
			case 'G':
				matches = true;
				break;
			}
			break;
		default:
			throw new Exception(String.format("Don't recognize IUPAC code \"%s\"", matchNuc));
		}
		
		return matches;
	}
	

}
