package edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers;

import java.util.*;

import org.biojava.bio.seq.StrandedFeature;

import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPos;
import edu.usc.epigenome.genomeLibs.Counters.SymbolCounterStratified;
import edu.usc.epigenome.genomeLibs.ReadPos.StreamHandlers.RPHandlerSymbolCountsStratifyByCycle;

public class APHandlerWindowStatsCpGConcordance extends APHandlerWindowStats {

	public APHandlerWindowStatsCpGConcordance(int inWindSize) {
		super(inWindSize);
	}

	@Override
	public void finish() {
		super.finish();
	}

	@Override
	public void init() {
		super.init();
	}


	@Override
	public boolean streamWindow(AlignmentPos[] priorAps,
			AlignmentPos currentAp, AlignmentPos[] nextAps,
			Queue<AlignmentPos> apWind) {

		
		SymbolCounterStratified windCounts = new SymbolCounterStratified();
		// Add our own counts
		windCounts.addCounts(currentAp.getSnpCounterStratifiedByCycle(true));
		// And out neighbors
		for(AlignmentPos ap : apWind)
		{
			windCounts.addCounts( ap.getSnpCounterStratifiedByCycle(true) );
		}
		
		double lowConversion = windCounts.getConvertedFrac(RPHandlerSymbolCountsStratifyByCycle.getLowStratString());
		double highConversion = windCounts.getConvertedFrac(RPHandlerSymbolCountsStratifyByCycle.getMediumStratString());
		
//		
//		Queue<AlignmentPos> apWindCopy = apWind;
//			
//			System.err.println(
//				currentAp + 
//				"\t" + 
//				AlignmentPos.getRefTokens(priorAps) + 
//				"/" + 
//				currentAp.getRefToken() + 
//				"/" +
//				AlignmentPos.getRefTokens(nextAps) + 
//				"\t" + 
//				AlignmentPos.getRefTokens(apWindCopy) +
//				"\t" + 
//				AlignmentPos.getPosString(apWindCopy) +
//				"\t" +
//				String.format("%.2f/%.2f", lowConversion, highConversion)
////				"\t" + 
////				windCounts.oneLineOutput()
//				);
//			
			
//			if ( (lowConversion != Double.NaN) && (highConversion != Double.NaN) )
			if ( (lowConversion >= 0.0) && (highConversion >= 0.0) )
			{
				double diff = highConversion - lowConversion;
				if ((diff!=0.0) && (diff > -1.0) && (diff < 1.0))
				{
//				System.out.println(String.format("%.2f,%.2f", lowConversion, highConversion));
				System.out.println(String.format("%.2f", diff));
				}
			}
		
		return true;
	}

}
