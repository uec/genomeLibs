package edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers;

import com.googlecode.charts4j.Color;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;

public class CpgNonconversionSummarizer extends CpgSummarizer {

	
	
	
	/**
	 * 
	 */
	public CpgNonconversionSummarizer() {
		super();
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param inQuerier
	 */
	public CpgNonconversionSummarizer(MethylDbQuerier inQuerier) {
		super(inQuerier);
		// TODO Auto-generated constructor stub
	}

	@Override
	public void init(MethylDbQuerier inQuerier) {
		// TODO Auto-generated method stub
		super.init(inQuerier);
		
		this.setColor(Color.ORANGE);
		this.setNaturalMax(10.0);
		this.setNaturalMin(0.0);
		
		this.setLabel("5' nonconversion frac");
		this.setDesc("");
	}

	@Override
	public void streamCpg(Cpg cpg) {
		// TODO Auto-generated method stub
		super.streamCpg(cpg);

		// Impose the filters just in case.
		if (cpg.passesOppositeAFilterDefault())
		{
			if (cpg.totalReadsOpposite >= 5)
			{
				double frac = (double)cpg.cReadsNonconversionFilt/((double)cpg.cReadsNonconversionFilt+(double)cpg.cReads+(double)cpg.tReads);
				this.streamValue( frac * 100.0, cpg.getCpgWeight()); // Make it percent
			}
		}
	}

	@Override
	public void removeCpg(Cpg cpg) {
		// Follow the same rules as the streamer
		if (cpg.passesOppositeAFilterDefault())
		{
			if (cpg.totalReadsOpposite >= 5)
			{
				double val = 100.0 * (double)cpg.cReadsNonconversionFilt/((double)cpg.cReadsNonconversionFilt+(double)cpg.cReads+(double)cpg.tReads);
				//System.err.println("Removing cpg with weight: " + cpg.getCpgWeight());
				this.removeValue(val, cpg.getCpgWeight());
			}
		}
	}


}
