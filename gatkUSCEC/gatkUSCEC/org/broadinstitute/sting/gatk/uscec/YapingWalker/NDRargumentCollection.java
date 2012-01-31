package org.broadinstitute.sting.gatk.uscec.YapingWalker;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteArgumentCollection;

public class NDRargumentCollection extends BisulfiteArgumentCollection {

	public NDRargumentCollection() {
		// TODO Auto-generated constructor stub
	}
	
	@Argument(fullName = "nucleosome_position_window", shortName = "npw", doc = "define the basic nucleosome depletion window size(bp)", required = false)
    public int nucPosWindow = 147;
	
	@Argument(fullName = "nucleosome_linker_window", shortName = "nlw", doc = "define the basic nucleosome linker window size(bp)", required = false)
    public int nucLinkerWindow = 40;
	
	@Argument(fullName = "minimum_number_gch_in_window_has_methy_value", shortName = "mgn", doc = "minimum number of gch in window has methy value", required = false)
    public int minGchNum = 3;
	
	@Argument(fullName = "minimum_CT_depth_for_gch_in_window", shortName = "mcd", doc = "minimum CT reads depth for GCH inside window", required = false)
    public int minCTDepth = 3;
	
	@Argument(fullName = "minimum_number_gch_in_linker_window_has_methy_value", shortName = "mgnlw", doc = "minimum number of gch in linker window has methy value", required = false)
    public int minGchNumLinkerWindow = 2;
	
	@Argument(fullName = "minimum_CT_depth_for_gch_in_linker_window", shortName = "mcdlw", doc = "minimum CT reads depth for GCH insidelinker  window", required = false)
    public int minCTDepthLinkerWindow = 3;
	
	@Argument(fullName = "wig_output", shortName = "wo", doc = "wig File to which variants should be written", required = true)
    public String wigFile = null;
	
	@Argument(fullName = "minimum_gch_methy_for_ndr", shortName = "ndrThreshold", doc = "minimum GCH methylation value criteria to be NDR region", required = false)
    public double ndrThreshold = 0.4;
	
	@Argument(fullName = "minimum_gch_methy_diff_for_ndr", shortName = "ndrDiffThreshold", doc = "minimum GCH methylation value differences with adjacent window to be identified as NDR region", required = false)
    public double ndrDiffThreshold = 0.4;
	
	@Argument(fullName = "enable_stat_test_for_ndr", shortName = "statTest", doc = "enable KS test rather than hard threshold to detect NDR region", required = false)
    public boolean statTest = false;
	
	@Argument(fullName = "sig_threshold_for_test", shortName = "sigValue", doc = "significance threshold to detect NDR region", required = false)
    public double sigValue = 0.01;
	
	
	public NDRargumentCollection clone() {
		NDRargumentCollection nac = new NDRargumentCollection();
		nac.nucPosWindow = nucPosWindow;
		nac.nucLinkerWindow = nucLinkerWindow;
		nac.minGchNum = minGchNum;
		nac.minCTDepth = minCTDepth;
		nac.minGchNumLinkerWindow = minGchNumLinkerWindow;
		nac.minCTDepthLinkerWindow = minCTDepthLinkerWindow;
		nac.wigFile = wigFile;
		nac.ndrThreshold = ndrThreshold;
		nac.ndrDiffThreshold = ndrDiffThreshold;
		nac.statTest = statTest;
		nac.sigValue = sigValue;
		
		return nac;
	}

}
