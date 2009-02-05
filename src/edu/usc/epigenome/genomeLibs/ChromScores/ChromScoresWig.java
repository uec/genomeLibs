package edu.usc.epigenome.genomeLibs.ChromScores;


public class ChromScoresWig extends ChromScoresArray {

	public ChromScoresWig(String fn, String genome)
	throws Exception
	{
		super(genome);
		this.populateFromWig(fn);
	}
	
}
