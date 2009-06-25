package edu.usc.epigenome.scripts;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.PileupToTemplate;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosOptions;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerFeatAligner;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer;

public class PileupToAlignedByFeats extends PileupToTemplate {

	@Option(name="-featGtf",usage="feature Gtf")
	protected String featGtf = null;
	@Option(name="-doCensoring",usage="Use feature censoring")
	protected boolean doCensoring = false;
    @Option(name="-featWindSize",usage="window size, for feature alignment (default 1000)")
    protected int featWindSize = -1;
 
	/**
	 * @param args
	 */
    public static void main(String[] args)
    throws Exception
    {
    	new PileupToAlignedByFeats().doMain(args);
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
		apStreamer.add(new APHandlerFeatAligner(featGtf, featWindSize, doCensoring));
	}


	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.PileupToTemplate#defaultApOptions()
	 */
	@Override
	protected AlignmentPosOptions defaultApOptions() {
		AlignmentPosOptions apos = super.defaultApOptions();
		apos.onlyFirstCycle = true;
		return apos;
	}

	
	
	
}
