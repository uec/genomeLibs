package edu.usc.epigenome.genomeLibs;

import org.biojava.bio.seq.DNATools;

public class APFilterCpgs extends AlignmentPosStreamFilter {

	public APFilterCpgs() {
		// TODO Auto-generated constructor stub
	}

	@Override
	public boolean elementPasses(AlignmentPos[] priorAps,
			AlignmentPos currentAp, AlignmentPos[] nextAps) {

		if ((nextAps.length < 1) || (priorAps.length < 1))
		{
			System.err.println("APFilterCpgs must have at least bp of sequence context on either side");
		//	(new Exception()).printStackTrace();
			System.exit(0);
		}
		
		return ((currentAp.getRef().equals(DNATools.c())) &&
				(nextAps[0].getRef().equals(DNATools.g())));
	}
}
