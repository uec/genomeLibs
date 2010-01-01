package edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers;

import com.googlecode.charts4j.Color;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;

public class CpgDeaminationSummarizer extends CpgSummarizer {

	
	
	
	/**
	 * 
	 */
	public CpgDeaminationSummarizer() {
		super();
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param inQuerier
	 */
	public CpgDeaminationSummarizer(MethylDbQuerier inQuerier) {
		super(inQuerier);
		// TODO Auto-generated constructor stub
	}

	@Override
	protected void init(MethylDbQuerier inQuerier) {
		// TODO Auto-generated method stub
		super.init(inQuerier);
		
		this.setColor(Color.ORANGE);
		this.setNaturalMax(2.0);
		this.setNaturalMin(0.0);
		
		this.setLabel("Genetic CtoT");
		this.setDesc("Greater than 20% A on opposite strand (at least 2)");
	}

	@Override
	public void streamCpg(Cpg cpg) {
		// TODO Auto-generated method stub
		super.streamCpg(cpg);
		
		// Make it a percent
		this.streamValue( (cpg.passesOppositeAFilterDefault()) ? 0.0 : 100.0);
	}


}
