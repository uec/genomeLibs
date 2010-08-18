package edu.usc.epigenome.genomeLibs.ChromScores;

import java.util.*;

public class ChromScoresMap extends ChromScoresFast {

	public ChromScoresMap(String genome) {
		super(genome);
	}

	
	protected Object newChromArray(int chr_len)
	{
		return new TreeMap<Integer,Number>();
	}
	

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.ChromScores.ChromScoresFast#setScore(java.lang.Object, int, java.lang.Number)
	 */
	@Override
	protected void setScore(Object array, int pos, Number score) {
		// TODO Auto-generated method stub
		Map<Integer,Number> m = (Map<Integer,Number>)array;
		Integer key = new Integer(pos);
		m.put(key, score);
	}


	protected Object addScoreToArray(Object array, int pos, Number score)
	{
		Map<Integer,Number> m = (Map<Integer,Number>)array;
		Integer key = new Integer(pos);
		Number oldVal = m.get(key);
		Number newVal = score;
		if (oldVal != null) newVal = new Double(oldVal.doubleValue() + newVal.doubleValue());
		m.put(key, newVal);
		return m;
	}
	
	protected Number getArrayScore(Object array, int pos)
	{
		Map<Integer,Number> m = (Map<Integer,Number>)array;
		Number out = (Number)m.get(new Integer(pos));
		
		// Should we auto this to 0?
		
		
		return out;
	}

		protected double[] getAllScores(Object array)
		{	
			return (double[])array;
		}

	 protected int minPos(Object array)
	 {
		 TreeMap<Integer,Number> m = (TreeMap<Integer,Number>)array;
		 int min_pos = ((Integer)m.firstKey()).intValue();
		 return min_pos;
	 }

	
	 protected int maxPos(Object array)
	 {
		 TreeMap<Integer,Number> m = (TreeMap<Integer,Number>)array;
		 int max_pos = ((Integer)m.lastKey()).intValue();
		 return max_pos;
	 }

	
	
}
