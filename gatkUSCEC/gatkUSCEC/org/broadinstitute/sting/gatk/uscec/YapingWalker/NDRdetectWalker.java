package org.broadinstitute.sting.gatk.uscec.YapingWalker;

import java.io.File;
import java.util.Iterator;
import java.util.LinkedList;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.DownsampleType;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;

import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisSNPUtils;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteDiploidSNPGenotypePriors;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteVariantCallContext;
import org.broadinstitute.sting.gatk.walkers.BAQMode;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.Downsample;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.baq.BAQ;
import org.broadinstitute.sting.utils.pileup.AbstractReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

@BAQMode(QualityMode = BAQ.QualityMode.OVERWRITE_QUALS, ApplicationTime = BAQ.ApplicationTime.ON_INPUT)
@Reference(window=@Window(start=-200,stop=200))
@By(DataSource.REFERENCE)
@Downsample(by=DownsampleType.NONE)
public class NDRdetectWalker extends LocusWalker<NDRCallContext,LinkedList<NDRCallContext>> implements
		TreeReducible<LinkedList<NDRCallContext>> {

	@ArgumentCollection private static NDRargumentCollection NAC = new NDRargumentCollection();
	
	
    protected WigWriterImp writer = null;
	
	//private NDRdetectionEngine NDRD_engine = null;
	
	private LinkedList<NDRCallContext> window = null;
	
	private GenotypePriors genotypePriors;
	
	private ContextCondition summary = null;
	
	public NDRdetectWalker() {
		// TODO Auto-generated constructor stub
		
	}

	 public static class ContextCondition {
		 	long nWindowsVisited = 0;

	        /** The number of bases that were potentially callable -- i.e., those not at excessive coverage or masked with N */
	     //   long nBasesCallable = 0;

	        /** The number of bases called confidently (according to user threshold), either ref or other */
	      //  long nBasesCalledConfidently = 0;

	        /** The number of bases for which calls were emitted */
	   //     long nCallsMade = 0;
	        
	        /** The average sequence depth inside window */
	        double sumGchCTDepthWind = 0;
	     
	        /** The number of Gch bases called confidently (according to user threshold), either ref or other */
	        long sumGchNumInWindCalledConfidently = 0;
	        
	        /** The number of Hcg bases called confidently (according to user threshold), either ref or other */
	  //      long nHcgBasesCalledConfidently = 0;
	        
	        /** The number of Wcg bases called confidently (according to user threshold), either ref or other */
	  //      long nWcgBasesCalledConfidently = 0;
	        
	        /** The sum of methylation value of Gch bases called confidently (according to user threshold), either ref or other */
	     //   double sumMethyGchBasesCalledConfidently = 0;
	        
	        /** The sum of methylation value of Hcg bases called confidently (according to user threshold), either ref or other */
	 //       double sumMethyHcgBasesCalledConfidently = 0;
	        
	        /** The sum of methylation value of Wcg bases called confidently (according to user threshold), either ref or other */
	 //       double sumMethyWcgBasesCalledConfidently = 0;
	        
	        /** The number of windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        long nWindowsCalledConfidently = 0;
	        
	        /** The number of NDR windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        long nNDRWindowsCalledConfidently = 0;
	        
	        /** The sum of methylation value of windows called confidently (according to user threshold), contained enough confidantly GCH and enough seq depth */
	        double sumMethyWindowsCalledConfidently = 0;
	        
	        double percentConfidantWindowOfAll() { return (double)nWindowsCalledConfidently/nWindowsVisited;}
	        double percentMethyConfidantWindowOfAll() { return (double)sumMethyWindowsCalledConfidently/(nWindowsCalledConfidently);}
	        double percentGchCTDepthConfidantWindowOfAll() { return (double)sumGchCTDepthWind/(nWindowsCalledConfidently);}
	        double percentGchNumInConfidantWindowOfAll() { return (double)sumGchNumInWindCalledConfidently/(nWindowsCalledConfidently);}
	        
	 }
	 
	 public void initialize() {
	//	 NDRD_engine = new NDRdetectionEngine(getToolkit(), NAC, logger);
		 genotypePriors = new BisulfiteDiploidSNPGenotypePriors();
		 File fn = new File(NAC.bedFile);
		 writer = new WigWriterImp(fn);
		 summary = new ContextCondition();
	 }
	
	@Override
	public NDRCallContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
		// TODO Auto-generated method stub
		String cytosinePattern = "GCH-2";
		double methyStatus = 0.5 ;
		BisSNPUtils it = new BisSNPUtils(NAC);
		AlignmentContext stratifiedContexts = it.getFilteredAndStratifiedContexts(NAC, ref, context);
		
		NDRCallContext ncc = new NDRCallContext(context, ref.getLocus());
		
		if(stratifiedContexts == null){
			ncc.setCytosinePatternFlag(false);
			return ncc;
		}
			
		
		if(it.checkCytosineStatus(cytosinePattern, stratifiedContexts.getBasePileup(), tracker, ref, (BisulfiteDiploidSNPGenotypePriors) genotypePriors, NAC, methyStatus)){
			//System.out.println(ref.getLocus().getStart());
			 ncc.setCytosinePatternFlag(true);
			 return ncc;
		}	
		else{
			//make fake context in chrM, so we know it immediately it is not Gch
			ncc.setCytosinePatternFlag(false);
			return ncc;
		}
		//return NDRD_engine.calculateNDRscore(tracker, ref, context);
			
	}

	@Override
	public LinkedList<NDRCallContext> reduce(NDRCallContext value, LinkedList<NDRCallContext> window) {
		// TODO Auto-generated method stub
		if ( value == null )
            return window;
		if(!window.isEmpty()){
			if(!value.getLoc().onSameContig(window.getFirst().getLoc())){
				window.clear();
				window.addLast(value);
				return window;
			}
			else if(value.getLoc().discontinuousP(window.getLast().getLoc())){
				window.clear();
				window.addLast(value);
				return window;
			}
		}
		
		
		if(window.size() <= NAC.nucPosWindow - 2){
			window.addLast(value);
			return window;
		}	
		else if(window.size() == NAC.nucPosWindow - 1){
			window.addLast(value);
		}
		else{
			window.removeFirst();
			window.addLast(value);
		}
		summary.nWindowsVisited++;
		double averageGchMethy = getValueFromWindow(window);
		//GenomeLoc locMidInWind = window.get(NAC.nucPosWindow/2).getLocation();
		if(!Double.isNaN(averageGchMethy)){
			summary.nWindowsCalledConfidently++;
			summary.sumMethyWindowsCalledConfidently += averageGchMethy;
			window.getLast().setGchMethyInWindow(averageGchMethy);
			if(averageGchMethy >= NAC.ndrThreshold){
				summary.nNDRWindowsCalledConfidently++;
				writer.add(window.getLast().getLoc(), averageGchMethy);
			}
			
			System.out.println(window.getFirst().getLoc().getStart() + "\t" + averageGchMethy);
		}
		
		return window;
	}

	@Override
	public LinkedList<NDRCallContext> reduceInit() {
		// TODO Auto-generated method stub
		window = new LinkedList<NDRCallContext>();
		return window;
	}
	
	@Override
	public LinkedList<NDRCallContext> treeReduce(LinkedList<NDRCallContext> lhWinds,
			LinkedList<NDRCallContext> rhWinds) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void onTraversalDone(LinkedList<NDRCallContext> sum) {
		logger.info(String.format("Visited windows                                %d", summary.nWindowsVisited));
		logger.info(String.format("Confidantly called windows                                %d", summary.nWindowsCalledConfidently));
		logger.info(String.format("Confidantly called NDR windows                                %d", summary.nNDRWindowsCalledConfidently));
		logger.info(String.format("Percentage of confidant called windows                                %.2f", summary.percentConfidantWindowOfAll()));
		logger.info(String.format("Average methylation in confidant called windows                                %.2f", summary.percentMethyConfidantWindowOfAll()));
		logger.info(String.format("Average GCH seq depth in confidant called windows                                %.2f", summary.percentGchCTDepthConfidantWindowOfAll()));
		logger.info(String.format("Average GCH number in confidant called windows                                %.2f", summary.percentGchNumInConfidantWindowOfAll()));
		writer.close();
	}

	public double getValueFromWindow(LinkedList<NDRCallContext> window){
		double averageGchMethy = Double.NaN;
		double sumGchMethy = 0;
		int numValidGch = 0;
		int sumGchCTDepth = 0;
		Iterator<NDRCallContext> itContext = window.iterator();
		while(itContext.hasNext()){
			NDRCallContext tmpContext = itContext.next();
			if(!tmpContext.getCytosinePatternFlag())
				continue;
			int numC = 0;
			int numT = 0;
			for( PileupElement p : tmpContext.getRealContext().getBasePileup()){
				boolean negStrand = p.getRead().getReadNegativeStrandFlag();
				if(((GATKSAMRecord)p.getRead()).isGoodBase(p.getOffset())){
					if(negStrand){
						if(p.getBase()==BaseUtils.G){
							numC++;
						}
						else if(p.getBase()==BaseUtils.A){
							numT++;
						}
						else{
							
						}
						
					}
					else{
						if(p.getBase()==BaseUtils.C){
							numC++;
						}
						else if(p.getBase()==BaseUtils.T){
							numT++;
						}
						else{
							
						}
					}
				}
				//if ( BisSNPUtils.usableBase(p, true) ){
					//sumGchSeqDepth++;
				//	if(p.getBase() == BaseUtils.C)
				//		numC++;
				//	if(p.getBase() == BaseUtils.T)
				//		numT++;
		        //}
			}
			if((numC + numT) >= NAC.minCTDepth){
				numValidGch++;
				sumGchCTDepth += (numC + numT);
				sumGchMethy += (double)numC/(double)(numC + numT);
			//	System.out.println("loc: " + tmpContext.getLoc().getStart() + "\tGchMethy: " + (double)numC/(double)(numC + numT) + "\tnumC: " + numC + "\tnumT: " + numT);
			}
		}
		
		if(numValidGch >= NAC.minGchNum){
			summary.sumGchNumInWindCalledConfidently += numValidGch;
			summary.sumGchCTDepthWind += sumGchCTDepth/(double)numValidGch;
			averageGchMethy = sumGchMethy/(double)numValidGch;
		}
			
		//System.out.println("numValidGch: " + numValidGch + "\tsumGchMethy: " + sumGchMethy + "\taverageGchMethy: " + averageGchMethy);
		
		return averageGchMethy;
	}

}
