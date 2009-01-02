package edu.usc.epigenome.genomeLibs;

import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.symbol.Symbol;

public class AlignmentPosSnpsBisulfiteConverted extends AlignmentPosSnps {

	
	public AlignmentPosSnpsBisulfiteConverted() {
		super();
	}

	public AlignmentPosSnpsBisulfiteConverted(char inRef, String inChr,
			int inPos, AlignmentPosOptions inApOptions) {
		super(inRef, inChr, inPos, inApOptions);
	}

	public AlignmentPosSnpsBisulfiteConverted(Symbol inRef, String inChr,
			int inPos, AlignmentPosOptions inApOptions) {
		super(inRef, inChr, inPos, inApOptions);
	}

	public AlignmentPosSnpsBisulfiteConverted(AlignmentPos ap) {
		super(ap);
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos#toGff(java.lang.String, boolean, int)
	 */
	@Override
	public SimpleGFFRecord toGff(boolean ref_fw_strand) {
		
		SimpleGFFRecord rec = super.toGff(ref_fw_strand);
		GFFUtils.add_gffrecord_map_entry(rec, "hiya", "1");
		return rec;
	}
	

}
