package org.broadinstitute.sting.gatk.examples;

import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDatumOptimized;
import org.broadinstitute.sting.utils.BaseUtils;

import org.broadinstitute.sting.gatk.examples.BisulfiteSnpUtil;

public class BisulfiteRecalDatumOptimized extends RecalDatumOptimized {

	public BisulfiteRecalDatumOptimized() {
		// TODO Auto-generated constructor stub
	}

	public BisulfiteRecalDatumOptimized(RecalDatumOptimized copy) {
		super(copy);
		// TODO Auto-generated constructor stub
	}

	public BisulfiteRecalDatumOptimized(long numObservations, long numMismatches) {
		super(numObservations, numMismatches);
		// TODO Auto-generated constructor stub
	}
	
	public synchronized final void incrementBaseCounts( final byte curBase, final byte refBase, final byte baseQual )  {
    	
		
            //if( readBaseIndex != refBaseIndex && !(BisulfiteSnpUtil.isCytosine(refBase,false) && BisulfiteSnpUtil.isCytosine(readBase,true))) {
                //counter.novelCountsMM++;
    		if( BaseUtils.simpleBaseToBaseIndex(curBase) != BaseUtils.simpleBaseToBaseIndex(refBase) ){
    			if((BisulfiteSnpUtil.isCytosine(refBase,false) && BisulfiteSnpUtil.isCytosine(curBase,true))){
    				if(baseQual <= 5){
    					increment( 1, 1 );
    				}
    				else{
    					increment( 0, 0 );
    				}
    			}
    			else{
    				increment( 1, 1 );
    			}
            }
    		else{
    			
    			increment( 1, 0 );
    		}


        increment( 1, (BisulfiteSnpUtil.isCytosine(refBase,false) && BisulfiteSnpUtil.isCytosine(curBase,true)) ? 0 : (BaseUtils.simpleBaseToBaseIndex(curBase) == BaseUtils.simpleBaseToBaseIndex(refBase) ? 0 : 1) ); // increment takes num observations, then num mismatches
    }
}
