package edu.usc.epigenome.genomeLibs.MethylDb;

import com.googlecode.charts4j.Color;

public class CpgCoverageSummarizer extends CpgSummarizer {

	
	
	
	/**
	 * 
	 */
	public CpgCoverageSummarizer() {
		super();
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param inQuerier
	 */
	public CpgCoverageSummarizer(MethylDbQuerier inQuerier) {
		super(inQuerier);
		// TODO Auto-generated constructor stub
	}

	@Override
	protected void init(MethylDbQuerier inQuerier) {
		// TODO Auto-generated method stub
		super.init(inQuerier);
		
		this.setColor(Color.GREEN);
		this.setNaturalMax(20);
		this.setNaturalMin(0.0);
		
		this.setLabel("Seq Depth");
		this.setDesc("Total reads on both strands");
	}

	@Override
	public void streamCpg(Cpg cpg) {
		// TODO Auto-generated method stub
		super.streamCpg(cpg);
		
		// Impose the filters just in case.
		
		
		double val = (double)(cpg.totalReads + cpg.totalReadsOpposite);
		this.streamValue(val);
	}


}
