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
		Iterator<AlignmentPos> it = new MaqPileupRealignmentIterator(fn, apos);
		
		while (it.hasNext())
		{
			AlignmentPos ap = it.next();
			System.err.println(ap.toString());
		}
	}

}
