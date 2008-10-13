package edu.usc.epigenome.genomeLibs;



public class FeatureCountParams {

	public int pre_gran, core_gran, post_gran;
	public int pre_size, post_size;

	
	public FeatureCountParams(int in_pre_size, int in_post_size,
			int in_pre_gran, int in_core_gran, int in_post_gran)
	{
		pre_size = in_pre_size;
		post_size = in_post_size;
		pre_gran = in_pre_gran;
		core_gran = in_core_gran;
		post_gran = in_post_gran;
	}
	
	public String toString()
	{
		return "pre=" + pre_size + ", post=" + post_size;
	}
	
}
