package edu.usc.epigenome.genomeLibs;

public class SymbolCounterStratifiedByCycle extends SymbolCounterStratified {

	private static final long serialVersionUID = 8373256943666483769L;
	
	protected int lowCycleStart = 1;
	protected int lowCycleEnd = 5;
	protected int mediumCycleStart = 16;
	protected int mediumCycleEnd = 20;
	protected int highCycleStart = 31;
	protected int highCycleEnd = 35;
	
	public SymbolCounterStratifiedByCycle() {
	}

	/**
	 * @param lowCycleStart
	 * @param lowCycleEnd
	 * @param highCycleStart
	 * @param highCycleEnd
	 */
	public SymbolCounterStratifiedByCycle(int lowCycleStart,
			int lowCycleEnd, int mediumCycleStart, int mediumCycleEnd, 
			int highCycleStart, int highCycleEnd) {
		super();
		this.lowCycleStart = lowCycleStart;
		this.lowCycleEnd = lowCycleEnd;
		this.mediumCycleStart = lowCycleStart;
		this.mediumCycleEnd = lowCycleEnd;
		this.highCycleStart = highCycleStart;
		this.highCycleEnd = highCycleEnd;
	}	
	

	/*** Getters ***/
	
	public int getLowCycleStart() {
		return lowCycleStart;
	}

	public int getLowCycleEnd() {
		return lowCycleEnd;
	}

	public int getMediumCycleStart() {
		return mediumCycleStart;
	}

	public int getMediumCycleEnd() {
		return mediumCycleEnd;
	}

	public int getHighCycleStart() {
		return highCycleStart;
	}

	public int getHighCycleEnd() {
		return highCycleEnd;
	}

	
	
	
	/*** Static stratification info ***/

	static protected String getLowStratString()
	{
		return "low";
	}
	
	static protected String getMediumStratString()
	{
		return "medium";
	}

	static protected String getHighStratString()
	{
		return "high";
	}
		
}
