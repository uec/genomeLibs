package edu.usc.epigenome.genomeLibs.BisulfiteSnp;

import org.broad.tribble.util.variantcontext.Allele;
import org.broadinstitute.sting.gatk.walkers.genotyper.BiallelicGenotypeLikelihoods;

public class BisulfiteBiallelicGenotypeLikelihoods{
    private String sample;
    private double[] GLs;
    private Allele A, B, C;
    private int depth;

    /**
     * Create a new object for sample with given alleles and genotype likelihoods
     *
     * @param sample              sample name
     * @param A                   allele A
     * @param B                   allele B
     * @param C                   allele C
     * @param log10AALikelihoods  AA likelihoods
     * @param log10ABLikelihoods  AB likelihoods
     * @param log10BBLikelihoods  BB likelihoods
     * @param log10BCLikelihoods  BB likelihoods
     * @param depth               the read depth used in creating the likelihoods
     */
    public BisulfiteBiallelicGenotypeLikelihoods(String sample,
                                        Allele A,
                                        Allele B,
                                        Allele C,
                                        double log10AALikelihoods,
                                        double log10ABLikelihoods,
                                        double log10BBLikelihoods,
                                        double log10BCLikelihoods,
                                        int depth) {
        this.sample = sample;
        this.A = A;
        this.B = B;
        this.C = C;
        this.GLs = new double[]{log10AALikelihoods, log10ABLikelihoods, log10BBLikelihoods, log10BCLikelihoods};
        this.depth = depth;
    }

    public String getSample() {
        return sample;
    }

    public double getAALikelihoods() {
        return GLs[0];
    }

    public double getABLikelihoods() {
        return GLs[1];
    }

    public double getBBLikelihoods() {
        return GLs[2];
    }
    
    public double getBCLikelihoods() {
        return GLs[3];
    }

    public double[] getLikelihoods() {
        return GLs;
    }

    public Allele getAlleleA() {
        return A;
    }

    public Allele getAlleleB() {
        return B;
    }

    public Allele getAlleleC() {
        return C;
    }
    
    public int getDepth() {
        return depth;
    }
    
}
