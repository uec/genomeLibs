package edu.usc.epigenome.genomeLibs.MethylDb;

import com.googlecode.charts4j.Color;

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
		
		// Impose the filters just in case.
		
		if (cpg.totalReadsOpposite >= 5)
		{
			double fracA = (double)cpg.aReadsOpposite/(double)cpg.totalReadsOpposite;
			//if (fracA>0.0) System.err.println(fracA);
			this.streamValue( (fracA>0.2) ? 100.0 : 0.0); // Make it percent
		}
	}


}
