package org.broadinstitute.sting.gatk.uscec.YapingWalker;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteArgumentCollection;

public class NDRargumentCollection extends BisulfiteArgumentCollection {

	public NDRargumentCollection() {
		// TODO Auto-generated constructor stub
	}
	
	@Argument(fullName = "nucleosome_position_window", shortName = "npw", doc = "define the basic nucleosome depletion window size(bp)", required = false)
    public int nucPosWindow = 150;
	
	@Argument(fullName = "minimum_number_gch_in_window_has_methy_value", shortName = "mgn", doc = "minimum number of gch in window has methy value", required = false)
    public int minGchNum = 5;
	
	@Argument(fullName = "minimum_CT_depth_for_gch_in_window", shortName = "mcd", doc = "minimum CT reads depth for GCH inside window", required = false)
    public int minCTDepth = 3;
	
	@Argument(fullName = "wig_output", shortName = "wo", doc = "wig File to which variants should be written", required = true)
    public String wigFile = null;
	
	@Argument(fullName = "minimum_gch_methy_for_ndr", shortName = "ndrThreshold", doc = "minimum GCH methylation value criteria to be NDR region", required = false)
    public double ndrThreshold = 0.4;
	
	@Argument(fullName = "minimum_gch_methy_diff_for_ndr", shortName = "ndrDiffThreshold", doc = "minimum GCH methylation value differences with adjacent window to be identified as NDR region", required = false)
    public double ndrDiffThreshold = 0.4;
	
	@Argument(fullName = "enable_ks_test_for_ndr", shortName = "ksTest", doc = "enable KS test rather than hard threshold to detect NDR region", required = false)
    public boolean ksTest = false;
	
	
	public NDRargumentCollection clone() {
		NDRargumentCollection nac = new NDRargumentCollection();
		nac.nucPosWindow = nucPosWindow;
		nac.minGchNum = minGchNum;
		nac.minCTDepth = minCTDepth;
		nac.wigFile = wigFile;
		nac.ndrThreshold = ndrThreshold;
		nac.ndrDiffThreshold = ndrDiffThreshold;
		nac.ksTest = ksTest;
		
		return nac;
	}

}
