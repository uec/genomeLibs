package org.broadinstitute.sting.gatk.examples;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.recalibration.CountCovariatesWalker;
import org.broadinstitute.sting.gatk.walkers.recalibration.Covariate;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalDataManager;
import org.broadinstitute.sting.gatk.walkers.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.recalibration.CountCovariatesWalker.CountedData;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.collections.NestedHashMap;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;

import org.broadinstitute.sting.gatk.examples.BisulfiteRecalDatumOptimized;
import org.broadinstitute.sting.gatk.examples.BisulfiteSnpUtil;

public class BisulfiteCountCovariatesWalker extends CountCovariatesWalker {
 
	
	static CountCovariatesWalker s = null;

		public BisulfiteCountCovariatesWalker() {
			// TODO Auto-generated constructor stub
			s = new CountCovariatesWalker();

		}


	/**
	 * 
     * Update the mismatch / total_base counts for a given class of loci.
     *
     * @param counter The CountedData to be updated
     * @param context The AlignmentContext which holds the reads covered by this locus
     * @param refBase The reference base
     */
	protected static void updateMismatchCounts(CountedData counter, final AlignmentContext context, final byte refBase) {
        for( PileupElement p : context.getBasePileup() ) {
            final byte readBase = p.getBase();
            final byte readBaseQual = p.getQual();
            final int readBaseIndex = BaseUtils.simpleBaseToBaseIndex(readBase);
            final int refBaseIndex  = BaseUtils.simpleBaseToBaseIndex(refBase);
            try{
            	if( readBaseIndex != -1 && refBaseIndex != -1 ) {
                    //if( readBaseIndex != refBaseIndex && !(BisulfiteSnpUtil.isCytosine(refBase,false) && BisulfiteSnpUtil.isCytosine(readBase,true))) {
                        //counter.novelCountsMM++;
            		if( readBaseIndex != refBaseIndex ){
            			if((BisulfiteSnpUtil.isCytosine(refBase,false) && BisulfiteSnpUtil.isCytosine(readBase,true))){
            				if(readBaseQual <= 5){
            					increaseNovelCountsMM(counter,1);
            				}
            				else{
            					continue;
            				}
            			}
            			else{
            				increaseNovelCountsMM(counter,1);
            			}
                    }
                    increaseNovelCountsBases(counter,1);
            	}
            }
            catch (InstantiationException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
			} catch (IllegalAccessException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
			}
        }
    }

    /**
     * 
     * Major workhorse routine for this walker.
     * Loop through the list of requested covariates and pick out the value from the read, offset, and reference
     * Using the list of covariate values as a key, pick out the RecalDatum and increment,
     *   adding one to the number of observations and potentially one to the number of mismatches
     * Lots of things are passed as parameters to this method as a strategy for optimizing the covariate.getValue calls
     *   because pulling things out of the SAMRecord is an expensive operation.
     * @param counter Data structure which holds the counted bases
     * @param gatkRead The SAMRecord holding all the data for this read
     * @param offset The offset in the read for this locus
     * @param refBase The reference base at this locus
     */
	@Override
	protected void updateDataFromRead(CountedData counter, final GATKSAMRecord gatkRead, final int offset, final byte refBase) {
        Object[][] covars;
		try {
			covars = (Comparable[][]) gatkRead.getTemporaryAttribute(getCovarsAttribute());
			// Using the list of covariate values as a key, pick out the RecalDatum from the data HashMap
	        NestedHashMap data;
	        data = dataManager.data;
	        final Object[] key = covars[offset];
			BisulfiteRecalDatumOptimized datum = (BisulfiteRecalDatumOptimized) data.get( key );
	        if( datum == null ) { // key doesn't exist yet in the map so make a new bucket and add it
	            // initialized with zeros, will be incremented at end of method
	            datum = (BisulfiteRecalDatumOptimized)data.put( new BisulfiteRecalDatumOptimized(), true, (Object[])key );
	        }

	        // Need the bases to determine whether or not we have a mismatch
	        final byte base = gatkRead.getReadBases()[offset];
	        final long curMismatches = datum.getNumMismatches();
	        final byte baseQual =  gatkRead.getBaseQualities()[offset];

	        // Add one to the number of observations and potentially one to the number of mismatches
	        datum.incrementBaseCounts( base, refBase, baseQual );
	        increaseCountedBases(counter,1);
			increaseNovelCountsBases(counter,1);
			increaseNovelCountsMM(counter,(datum.getNumMismatches() - curMismatches));
	        
		} catch (IllegalArgumentException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IllegalAccessException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (InstantiationException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
        //counter.novelCountsBases++;
        
        //counter.novelCountsMM += datum.getNumMismatches() - curMismatches; // For sanity check to ensure novel mismatch rate vs dnsnp mismatch rate is reasonable
       
    }

	@Override
		protected void printMappingsSorted( final PrintStream recalTableStream, final int curPos, final Object[] key, final Map data) {
	        final ArrayList<Comparable> keyList = new ArrayList<Comparable>();
	        for( Object comp : data.keySet() ) {
	            keyList.add((Comparable) comp);
	        }

	        Collections.sort(keyList);
			//recalTableStream.println("haha" );
	        for( Comparable comp : keyList ) {
	            key[curPos] = comp;
	            final Object val = data.get(comp);
	            if( val instanceof BisulfiteRecalDatumOptimized ) { // We are at the end of the nested hash maps
	                // For each Covariate in the key
	                for( Object compToPrint : key ) {
	                    // Output the Covariate's value
	                    recalTableStream.print( compToPrint + "," );
	                }
	                // Output the RecalDatum entry
	                recalTableStream.println( ((BisulfiteRecalDatumOptimized)val).outputToCSV() );
	            } else { // Another layer in the nested hash map
	                printMappingsSorted( recalTableStream, curPos + 1, key, (Map) val );
	            }
	        }
	    }

		@Override
		protected void printMappings( final PrintStream recalTableStream, final int curPos, final Object[] key, final Map data) {
	        for( Object comp : data.keySet() ) {
	            key[curPos] = comp;
	            final Object val = data.get(comp);
	            if( val instanceof BisulfiteRecalDatumOptimized ) { // We are at the end of the nested hash maps
	                // For each Covariate in the key
	                for( Object compToPrint : key ) {
	                    // Output the Covariate's value
	                    recalTableStream.print( compToPrint + "," );
	                }
	                // Output the RecalDatum entry
	                recalTableStream.println( ((BisulfiteRecalDatumOptimized)val).outputToCSV() );
	            } else { // Another layer in the nested hash map
	                printMappings( recalTableStream, curPos + 1, key, (Map) val );
	            }
	        }
	    }
    
    private static void increaseCountedData(CountedData counter, int i, long increaseNum) throws InstantiationException, IllegalAccessException {
    	//CountedData s = CountedData.class.newInstance();     
    	Field f[] = counter.getClass().getDeclaredFields(); 
    	f[i].setAccessible(true); 
    	long num = f[i].getLong(counter);
    	num = num + increaseNum;
    	f[i].setLong(counter, num);
		//System.out.println("Here is your " + f[i].getName() + " = " + num );
    }
    
    private static void increaseNovelCountsMM(CountedData counter, long increaseNum) throws InstantiationException, IllegalAccessException {
    	increaseCountedData(counter, 7, increaseNum);
    }
    
    private static void increaseNovelCountsBases(CountedData counter, long increaseNum) throws InstantiationException, IllegalAccessException {
    	increaseCountedData(counter, 8, increaseNum);
    }
    
    private static void increaseCountedBases(CountedData counter, long increaseNum) throws InstantiationException, IllegalAccessException {
    	increaseCountedData(counter, 1, increaseNum);
    }
    
 /*   private static void increaseCountedSites(CountedData counter, long increaseNum) throws InstantiationException, IllegalAccessException {
    	increaseCountedData(counter, 0, increaseNum);
    }
    
    private static void increaseSkippedSites(CountedData counter, long increaseNum) throws InstantiationException, IllegalAccessException {
    	increaseCountedData(counter, 2, increaseNum);
    }
    
    private static void increaseSolidInsertedReferenceBases(CountedData counter, long increaseNum) throws InstantiationException, IllegalAccessException {
    	increaseCountedData(counter, 3, increaseNum);
    }
    
    private static void increaseOtherColorSpaceInconsistency(CountedData counter, long increaseNum) throws InstantiationException, IllegalAccessException {
    	increaseCountedData(counter, 4, increaseNum);
    }*/
    
    private static String getCovarsAttribute() throws IllegalArgumentException, IllegalAccessException, InstantiationException{
    	//CountCovariatesWalker s = CountCovariatesWalker.class.newInstance();     
    	Field f[] = s.getClass().getDeclaredFields(); 
    	f[2].setAccessible(true);     
	    return f[2].get(s).toString(); 
    }
    

 /*   private static String getSkipRecordAttribute() throws IllegalArgumentException, IllegalAccessException, InstantiationException{
    	CountCovariatesWalker s = CountCovariatesWalker.class.newInstance();     
    	Field f[] = s.getClass().getDeclaredFields(); 
    	f[0].setAccessible(true);     
	    return f[0].get(s).toString(); 
    }
    
    private static String getSeenAttribute() throws IllegalArgumentException, IllegalAccessException, InstantiationException{
    	CountCovariatesWalker s = CountCovariatesWalker.class.newInstance();     
    	Field f[] = s.getClass().getDeclaredFields(); 
    	f[1].setAccessible(true);     
	    return f[1].get(s).toString(); 
    }*/
    
 /*   private static RecalDataManager getDataManger() throws IllegalArgumentException, IllegalAccessException, InstantiationException{
    	//CountCovariatesWalker s = CountCovariatesWalker.class.newInstance();     
    	Field f[] = s.getClass().getDeclaredFields(); 
    	f[11].setAccessible(true);     
	    return (RecalDataManager) f[11].get(s); 
    }

    private static RecalibrationArgumentCollection getRac() throws IllegalArgumentException, IllegalAccessException, InstantiationException{
    	CountCovariatesWalker s = CountCovariatesWalker.class.newInstance();     
    	Field f[] = s.getClass().getDeclaredFields(); 
    	f[3].setAccessible(true);     
	    return (RecalibrationArgumentCollection) f[3].get(s); 
    }

    private static ArrayList<Covariate> getRequestedCovariates() throws IllegalArgumentException, IllegalAccessException, InstantiationException{
    	CountCovariatesWalker s = CountCovariatesWalker.class.newInstance();     
    	Field f[] = s.getClass().getDeclaredFields(); 
    	f[12].setAccessible(true);     
	    return (ArrayList<Covariate>) f[12].get(s); 
    }
    */
}
