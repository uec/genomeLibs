package edu.usc.epigenome.genomeLibs;

import java.util.Iterator;

public class ChromScoresIteratorAlignmentPos extends ChromScoresIterator {

	protected Iterator<AlignmentPos> f_ap_iter = null;
	protected String f_last_chr = null;
	protected int f_last_pos = -1;
	protected ChromScoresFast f_chrom_scores= null;
	String f_genome = null;
	
	public ChromScoresIteratorAlignmentPos(Iterator<AlignmentPos> ap_iter,
			String genome)
	{
		// Default method goes through chromosome by chromosome
		f_ap_iter = ap_iter;
		f_genome = genome;

		startNewChromScores();
	}
	
	
	@Override
	public boolean hasNext() {
		return f_ap_iter.hasNext();
	}

	@Override
	public ChromScoresFast next() {

		ChromScoresFast out = f_chrom_scores;
		
		// We go until we have no more APs or until we reach
		// one with a new chromosome.
		boolean done = false;
		while (!done && f_ap_iter.hasNext())
		{
			AlignmentPos ap = f_ap_iter.next();
			if (onNewChromScores(ap))
			{
				done = true;
				startNewChromScores();
			}

//			System.err.println("New ap: " + ap.f_chr + ": " + ap.f_pos);
			
			addAlignmentPos(ap);
			f_last_chr = ap.f_chr;
			f_last_pos = ap.f_pos;
		}
		
		return out;
	}
	
	protected void addAlignmentPos(AlignmentPos ap)
	{
		f_chrom_scores.addScore(ap.f_chr, ap.f_pos, ap.getTotalDepth());
	}
	
	
	// Returns true if the ap is the start of a new chromScores. 
	protected boolean onNewChromScores(AlignmentPos ap)
	{
//		System.err.println("Start new ChromScoresFast for chrom (" + ap.f_chr + "!=" + f_last_chr +")?");
		
		return (!((f_last_chr == null) || ap.f_chr.equals(f_last_chr)));
	}
	
	protected void startNewChromScores()
	{
		System.err.println("Starting new ChromScoresFast for chrom");
		f_chrom_scores = new ChromScoresArray(f_genome);
	}


}
