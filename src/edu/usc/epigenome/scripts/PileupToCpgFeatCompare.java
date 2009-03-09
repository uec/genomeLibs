package edu.usc.epigenome.scripts;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.PileupToCpgTemplate;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerCpgFeatAligner;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerCpgFeatComparer;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerCpgFilterNonSnpCpgs;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerCpgWindowAutocorr;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerCpgWindowAutocorrStranded;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer;

public class PileupToCpgFeatCompare extends PileupToCpgTemplate {

	@Option(name="-featGtf",usage="feature Gtf")
	protected String featGtf = null;
//	@Option(name="-featGtfs",usage="feature Gtfs")
//	protected String[] featGtfs = null;
    @Option(name="-featWindSize",usage="will compare to anything in this window (default 1000)")
    protected int featWindSize = -1;
 
	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new PileupToCpgFeatCompare().doMain(args);
    }

     
 	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.PileupToCpgTemplate#checkArgs()
	 */
	@Override
	public void checkArgs() throws Exception {
		super.checkArgs();

		if (featGtf == null) throw new CmdLineException("Must supply featGtf file");
		if (featWindSize == -1) throw new CmdLineException("Must supply featWindSize file");
	}


    /* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.PileupToCpgTemplate#addHandlers(edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer)
	 */
	@Override
	protected void addHandlers(AlignmentPosStreamer apStreamer) {
		apStreamer.add(new APHandlerCpgFeatComparer(featGtf, featWindSize));
	}

	
	
	
}
