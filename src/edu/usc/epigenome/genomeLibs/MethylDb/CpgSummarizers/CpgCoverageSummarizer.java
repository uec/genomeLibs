package edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers;

import com.googlecode.charts4j.Color;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;

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
	public void init(MethylDbQuerier inQuerier) {
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
		
		// FORWARD STRAND ONLY
		double val = (double)(cpg.totalReads);
		this.streamValue(val);
	}


}
