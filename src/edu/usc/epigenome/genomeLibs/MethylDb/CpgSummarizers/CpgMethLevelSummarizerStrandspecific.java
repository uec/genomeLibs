package edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.MethylDb.Cpg;
import edu.usc.epigenome.genomeLibs.MethylDb.MethylDbQuerier;

public class CpgMethLevelSummarizerStrandspecific extends
		CpgMethLevelSummarizer {

	protected boolean fwStrand = true; // If set to false, it's reverse strand
	
	public CpgMethLevelSummarizerStrandspecific(boolean inFwStrand) {
		fwStrand = inFwStrand;
	}

	public CpgMethLevelSummarizerStrandspecific(MethylDbQuerier inQuerier, boolean inFwStrand) {
		super(inQuerier);
		fwStrand = inFwStrand;
		
	}

	/**
	 * @return true if FW strand, false if REV strand
	 */
	public boolean isFwStrand() {
		return fwStrand;
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer#streamCpg(edu.usc.epigenome.genomeLibs.MethylDb.Cpg)
	 */
	@Override
	public void streamCpg(Cpg cpg) {

		StrandedFeature.Strand strand = cpg.getStrand();
		boolean use = ((isFwStrand() && (strand == StrandedFeature.POSITIVE)) ||
				(!isFwStrand() && (strand == StrandedFeature.NEGATIVE)));
				
		if (use)
		{
			super.streamCpg(cpg);
		}
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer#removeCpg(edu.usc.epigenome.genomeLibs.MethylDb.Cpg)
	 */
	@Override
	public void removeCpg(Cpg cpg) {
		StrandedFeature.Strand strand = cpg.getStrand();
		boolean use = ((isFwStrand() && (strand == StrandedFeature.POSITIVE)) ||
				(!isFwStrand() && (strand == StrandedFeature.NEGATIVE)));

		if (use)
		{
			super.removeCpg(cpg);
		}
	}

	
	
}
