package edu.usc.epigenome.genomeLibs;

public class MergeGffParams {

	public boolean f_combine_sources = false;
	public boolean f_combine_atts = false;
	public boolean f_combine_scores = false;
	
	public double f_min_score = Double.NEGATIVE_INFINITY; 
	public int f_min_features = Integer.MIN_VALUE; 
	
	public MergeGffParams()
	{
	}

	public MergeGffParams(boolean combine_sources, boolean combine_atts, 
			boolean combine_scores, double min_score, int min_features)
	{
		f_combine_sources = combine_sources;
		f_combine_atts = combine_atts;
		f_combine_scores = combine_scores;
		f_min_score = min_score; 
		f_min_features = min_features; 
	}

}
