package org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel;

import org.broad.tribble.util.variantcontext.Allele;
import org.broadinstitute.sting.gatk.walkers.genotyper.BiallelicGenotypeLikelihoods;

public class BisulfiteBiallelicGenotypeLikelihoods extends
		BiallelicGenotypeLikelihoods {

	private Integer[] CytosineStatus;
	
	public BisulfiteBiallelicGenotypeLikelihoods(String sample, Allele A,
			Allele B, double log10aaLikelihoods, double log10abLikelihoods,
			double log10bbLikelihoods, int depth) {
		super(sample, A, B, log10aaLikelihoods, log10abLikelihoods,
				log10bbLikelihoods, depth);
		// TODO Auto-generated constructor stub
	}
	
	public BisulfiteBiallelicGenotypeLikelihoods(String sample, Allele A,
			Allele B, double log10aaLikelihoods, double log10abLikelihoods,
			double log10bbLikelihoods, int depth, Integer[] CytosineStatus) {
		super(sample, A, B, log10aaLikelihoods, log10abLikelihoods,
				log10bbLikelihoods, depth);
		this.CytosineStatus = CytosineStatus;
		// TODO Auto-generated constructor stub
	}
	
	public Integer[] getCytosineStatus() {
        return this.CytosineStatus;
    }

}
