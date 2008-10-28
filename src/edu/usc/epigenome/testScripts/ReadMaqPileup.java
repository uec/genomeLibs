package edu.usc.epigenome.testScripts;

import java.util.Iterator;
import edu.usc.epigenome.genomeLibs.*;

public class ReadMaqPileup {

	private static final String C_USAGE = "Usage: ReadMaqPileup max_identical file1.pileup";
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
	throws Exception {

		if (args.length < 2)
		{
			System.err.println(C_USAGE);
			System.exit(0);
		}

		AlignmentPosOptions apos = new AlignmentPosOptions();
		String fn = null;
		try
		{
			apos.maxIdentical = Integer.parseInt(args[0]);
			fn = args[1];
		}
		catch (Exception e)
		{
			System.err.println(C_USAGE + "\n" + e);
			System.exit(0);
		}
			
		Iterator<AlignmentPos> ap_it = new AlignmentPosIteratorMaqPileup(fn, apos);
	
		while (ap_it.hasNext())
		{
			AlignmentPos ap = ap_it.next();
			System.out.println(ap.toString());
		}
	}
		
}
