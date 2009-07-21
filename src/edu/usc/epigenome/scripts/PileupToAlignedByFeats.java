package edu.usc.epigenome.scripts;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.PileupToTemplate;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosOptions;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerFeatAligner;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.APHandlerFeatAlignerEachfeat;
import edu.usc.epigenome.genomeLibs.AlignmentPos.StreamHandlers.AlignmentPosStreamHandler;
import edu.usc.epigenome.genomeLibs.AlignmentPos.Streamers.AlignmentPosStreamer;

public class PileupToAlignedByFeats extends PileupToTemplate {

	@Option(name="-featGtf",usage="feature Gtf")
	protected String featGtf = null;
	@Option(name="-doCensoring",usage="Use feature censoring")
	protected boolean doCensoring = false;
	@Option(name="-eachFeat",usage="Output a line for each feature (ouput can be huge, default false)")
	protected boolean eachFeat = false;
    @Option(name="-featWindSize",usage="window size, for feature alignment (default 1000)")
    protected int featWindSize = -1;
    @Option(name="-fragSize",usage="fragment size, for feature alignment (default 500)")
    protected int fragSize = 500;
    @Option(name="-downsamplingFactor",usage="Downsamples data matrix by this factor, for instance 10 will yield a 10x smaller matrix (default 1)")
    protected int downsamplingFactor = 1;
   
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
		
		if (!eachFeat && (downsamplingFactor != 1))
		{
			System.err.println("Can only use -downsamplingFactor with -eachFeat");
			System.exit(1);
		}
		
		AlignmentPosStreamHandler featAligner = (eachFeat) ? 
				new APHandlerFeatAlignerEachfeat(featGtf, featWindSize, doCensoring, fragSize,downsamplingFactor) : 
			new APHandlerFeatAligner(featGtf, featWindSize, doCensoring, fragSize);
		
		apStreamer.add(featAligner);
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
