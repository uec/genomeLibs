package edu.usc.epigenome.scripts;

import edu.usc.epigenome.genomeLibs.PileupToCpgTemplate;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerCpgFilterNonSnpCpgs;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerCpgWindowAutocorr;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerCpgWindowAutocorrStranded;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer;

public class PileupToCpgAutocorrelations extends PileupToCpgTemplate {

	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new PileupToCpgAutocorrelations().doMain(args);
    }

    /* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.PileupToCpgTemplate#addHandlers(edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer)
	 */
	@Override
	protected void addHandlers(AlignmentPosStreamer apStreamer) {
		
		apStreamer.add(new APHandlerCpgFilterNonSnpCpgs());
//		apStreamer.add(new APHandlerCpgWindowAutocorr(3000));
		apStreamer.add(new APHandlerCpgWindowAutocorrStranded(25000));

	}

	
	
	
}
