package edu.usc.epigenome.genomeLibs;

public class ChromScoresFastCounter extends ChromScoresArray {
//public class ChromScoresFastCounter extends ChromScoresMap {

	protected int f_ranges_added = 0;
	
	public ChromScoresFastCounter() {
		init();
	}

	
	public void addRange(String chr, int s, int e, int count)
	{
		addRange(chr, s, e, new Integer(count));
	}
	
//	public void addRange(String chr, int s, int e, Number score)
//	{
//		int chr_num = (new ChromFeatures()).chrom_from_public_str(chr);
//		
//		//System.err.println("Adding range " + chr + ":" + s + "-" + e + "\tcount=" + count);
//		for (int i = s; i <= e; i++)
//		{
//			addScore(chr, chr_num, i, score);
//		}
//		
//		f_ranges_added++;
//		if ((f_ranges_added % 1000) == 0)
//		{
//			System.err.println(f_ranges_added + "\tranges added");
//		}
//		
//	}
	
	protected Object addScoreToArray(Object array, int pos, Number score)
	{
		Number old_score = (Number)getArrayScore(array,pos);
		if (old_score == null) old_score = new Integer(0);
		super.addScoreToArray(array, pos, new Integer(old_score.intValue()+1));
		return array;
	}

	

	
	
}
