package edu.usc.epigenome.testScripts;

import java.io.PrintWriter;
import java.util.Iterator;
import edu.usc.epigenome.genomeLibs.*;

public class MaqPileupToWig {

	private static final String C_USAGE = "Usage: MaqPileupToWig file1.pileup name " + 
	" span step [view_limit_lower view_limit_upper]";
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception {

		if (args.length < 4)
		{
			System.err.println(C_USAGE);
			System.exit(0);
		}
		
		WigOptions wo = new WigOptions();

		String fn = args[0];
		String name = args[1];
		int span = Integer.parseInt(args[2]);
		int step = Integer.parseInt(args[3]);
		if (args.length > 4)
		{
			int view_lower = Integer.parseInt(args[4]);
			int view_upper = Integer.parseInt(args[5]);
			wo.f_view_limits = new int[] {view_lower, view_upper};
		}
		wo.f_span = span;
		wo.f_step = step;
		wo.f_format = 1;
		
		
		AlignmentPosOptions apos = new AlignmentPosOptions();
		Iterator<AlignmentPos> ap_it = new AlignmentPosIteratorMaqPileup(fn, apos);
		
		ChromScoresIteratorAlignmentPosFwRev cs_it = new ChromScoresIteratorAlignmentPosFwRev(ap_it, apos.f_genome);
		
		
		PrintWriter pw = new PrintWriter(System.out);
		while (cs_it.hasNext())
		{
			ChromScoresFast[] cs = cs_it.next();
			System.err.println("New chromScores iterator:" + cs);

			cs[0] = cs[0].smooth(wo.f_span, 36);
			wo.f_name = name + "fw";
			wo.makeFwStrand();
			cs[0].wigOutput(pw, wo);

			cs[1] = cs[1].smooth(wo.f_span, 36);
			wo.f_name = name + "rev";
			wo.makeRevStrand();
			cs[1].wigOutput(pw, wo);

		}
	}

}
