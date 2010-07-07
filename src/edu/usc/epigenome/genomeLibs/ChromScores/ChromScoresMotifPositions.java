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
		int lastPos = full_residues_arr.length - 1;
		Integer oneAsInt = new Integer(1);
		
		boolean doFw = (desiredStrand == StrandedFeature.POSITIVE) || (desiredStrand == StrandedFeature.UNKNOWN);
		boolean doRev = (desiredStrand == StrandedFeature.NEGATIVE) || (desiredStrand == StrandedFeature.UNKNOWN);
		
		char[] motifRevcompArr = MiscUtils.revCompNucStr(motif).toCharArray();
		
		for (int startPos = 0; startPos<=lastPos; startPos++)
		{
			if (doFw)
			{
				boolean match = true;
				for (int motifPos = 0; match && (motifPos<motifLen); motifPos++)
				{
					int refPos = motifPos + startPos;
					match = (full_residues_arr[refPos] == motifArr[motifPos]);
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
					match = (full_residues_arr[refPos] == motifRevcompArr[motifPos]);
				}
				if (match)
				{
					//System.err.printf("Found a (-) strand match to %s at pos %d\n", motif, start_coord+startPos+motifLen-1);
					this.addScore(target_chr, start_coord+startPos+motifLen-1, oneAsInt);
				}
			}
		}
	}
	

	

}
