package edu.usc.epigenome.scripts;

import edu.usc.epigenome.genomeLibs.PileupToCpgTemplate;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerCpgEmitter;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer;

public class PileupToCpgCsv extends PileupToCpgTemplate {

	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new PileupToCpgCsv().doMain(args);
    }

    /* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.PileupToCpgTemplate#addHandlers(edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer)
	 */
	@Override
	protected void addHandlers(AlignmentPosStreamer apStreamer) {
		
		System.err.println("Adding handler");
		apStreamer.add(new APHandlerCpgEmitter(this.cpgTrackFilename));

	}

	
	
	
}
