package edu.usc.epigenome.testScripts;

import java.util.Iterator;
import edu.usc.epigenome.genomeLibs.*;

public class ReadMaqPileup {

	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception {

		String fn = args[0];
		
		AlignmentPosOptions apos = new AlignmentPosOptions();
		Iterator<AlignmentPos> ap_it = new MaqPileupAlignmentPosIterator(fn, apos);
		
		ChromScoresIteratorAlignmentPos cs_it = new ChromScoresIteratorAlignmentPos(ap_it, apos.f_genome);
		
		while (cs_it.hasNext())
		{
			ChromScoresFast cs = cs_it.next();
			System.err.println("New chromScores iterator:" + cs);
			cs = cs.smooth(500, 36);
			cs.wigOutput("name_hiya", 0.0, "test.wig", true, 1, 500);
			
		}
	}

}
