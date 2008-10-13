package edu.usc.epigenome.genomeLibs;

public class ChromScoresWig extends ChromScoresArray {

	public ChromScoresWig(String fn)
	throws Exception
	{
		init();
		this.populateFromWig(fn);
	}

	public ChromScoresWig(String fn, String genome)
	throws Exception
	{
		init(genome);
		this.populateFromWig(fn);
	}

	
}
