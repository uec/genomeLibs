package edu.usc.epigenome.genomeLibs.AlignmentPos;

import java.util.*;

import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.biojava.bio.symbol.Symbol;

import edu.usc.epigenome.genomeLibs.GFFUtils;
import edu.usc.epigenome.genomeLibs.Counters.SymbolCounter;


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

	
	/****
	 * conversion counting (fw strand only)
	 */
	
	public double convertedFrac()
	{
		SymbolCounter snpCounts = this.getSnpCounter(true);
		double conv = snpCounts.getConvertedFrac();
		return conv;
	}
	
	

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPos#toGff(java.lang.String, boolean, int)
	 */
	@Override
	public SimpleGFFRecord toGff(boolean ref_fw_strand) {
		
		SimpleGFFRecord rec = super.toGff(ref_fw_strand);
		GFFUtils.add_gffrecord_map_entry(rec, "conversion", "" + this.convertedFrac());
		return rec;
	}
	
	
	static String getConvertedFracString(Collection<AlignmentPos> aps)
	{
		int len = aps.size();
		StringBuffer buf = new StringBuffer(len*5);
		Iterator<AlignmentPos> it = aps.iterator();
		while (it.hasNext())
		{
			AlignmentPos ap = it.next();
			AlignmentPosSnpsBisulfiteConverted apCast = (AlignmentPosSnpsBisulfiteConverted)ap;
			buf.append(String.format("%.2f,",apCast.convertedFrac()));
		}
		return buf.toString();
	}
		

}
