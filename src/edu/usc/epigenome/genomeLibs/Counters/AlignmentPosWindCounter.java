/**
 * 
 */
package edu.usc.epigenome.genomeLibs.Counters;

import java.io.PrintStream;
import java.util.Iterator;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
//import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosOptions;
import edu.usc.epigenome.genomeLibs.GenomicRange.GenomicRange;


/**
 * @author benb
 * 
 *
 */
public class AlignmentPosWindCounter extends GenomicRangeWindCounter  {


	private static final long serialVersionUID = -8312065178497620972L;

	/**
	 * 
	 */
	public AlignmentPosWindCounter(int inWindSize, boolean inStrandSpecific, String inGenome) {
		super(inWindSize, inStrandSpecific, inGenome);
	}

	
	public void increment(AlignmentPos ap) {
		
//		if (apOptions==null) apOptions = ap.getApOptions();
				
		if (strandSpecific)
		{
			super.increment(ap.getGenomicRangeStranded(StrandedFeature.POSITIVE), ap.getDepth(true));
			super.increment(ap.getGenomicRangeStranded(StrandedFeature.NEGATIVE), ap.getDepth(false));
		}
		else
		{
			super.increment(ap.getGenomicRangeStranded(StrandedFeature.UNKNOWN), ap.getTotalDepth());
		}
		
		//System.err.println("Counter keys: " + this.size() + " (total count " + this.getTotalCount() + ")");
	}

	
}
