package edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers;

import com.googlecode.charts4j.Color;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;

public class CpgMethLevelSummarizer extends CpgSummarizer {

	
	
	
	/**
	 * 
	 */
	public CpgMethLevelSummarizer() {
		super();
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param inQuerier
	 */
	public CpgMethLevelSummarizer(MethylDbQuerier inQuerier) {
		super(inQuerier);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void init(MethylDbQuerier inQuerier) {
		// TODO Auto-generated method stub
		super.init(inQuerier);
		
		this.setColor(Color.BLUE);
		this.setNaturalMax(1.0);
		this.setNaturalMin(0.0);
		
		this.setLabel("%Meth");
	}

	@Override
	public void streamCpg(Cpg cpg) {
		super.streamCpg(cpg);
		
		// Impose the filters just in case.
		if (cpg.passesOppositeAFilterDefault())
		{
			double val = cpg.fracMeth(true);
			this.streamValue(val);
		}
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgSummarizer#removeCpg(edu.usc.epigenome.genomeLibs.MethylDb.Cpg)
	 */
	@Override
	public void removeCpg(Cpg cpg) {
		// Follow the same rules as the streamer
		if (cpg.passesOppositeAFilterDefault())
		{
			double val = cpg.fracMeth(true);
			this.removeValue(val);
		}
	}
	
	
}
