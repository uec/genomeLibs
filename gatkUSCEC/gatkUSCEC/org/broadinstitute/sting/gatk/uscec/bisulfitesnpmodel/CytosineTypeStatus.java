package org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel;

import java.util.HashMap;

import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel;

public class CytosineTypeStatus {
	boolean isC = false;
	boolean isCpg = false;
	boolean isChh = false;
	boolean isChg = false;
	boolean isGch = false;
	boolean isGcg = false;
	boolean isHcg = false;
	double cytosineMethyLevel = 0;
	double cpgMethyLevel = 0;
	double chgMethyLevel = 0;
	double chhMethyLevel = 0;
	double gchMethyLevel = 0;
	double gcgMethyLevel = 0;
	double hcgMethyLevel = 0;
	
	BisulfiteArgumentCollection BAC = null;
	
	int maxCytosineLength = 3;
	
	HashMap<String, Double[]> cytosineListMap = null;
	
	public CytosineTypeStatus(BisulfiteArgumentCollection bac){ 
		
		BAC = bac;
		makeCytosineMap();
	}
	
	public void makeCytosineMap(){
		makeCytosineMap(BAC);
	}
	
	public void makeCytosineMap(BisulfiteArgumentCollection BAC){
		cytosineListMap = new HashMap<String, Double[]>();
		Double[] tmpDouble = new Double[4];//tmpDouble[0]: log10 likelihood in positive strand; tmpDouble[1]: log10 likelihood in negative strand; tmpDouble[2]: methylation level; tmpDouble[3]: is it this type of cytosine: 1-true or 0-false? 
		tmpDouble[0] = Double.NEGATIVE_INFINITY;
		tmpDouble[1] = Double.NEGATIVE_INFINITY;
		tmpDouble[2] = BAC.forceChh;
		tmpDouble[3] = 0.0;
		cytosineListMap.put("C-1".toUpperCase(), tmpDouble);
		tmpDouble = new Double[4];
		tmpDouble[0] = Double.NEGATIVE_INFINITY;
		tmpDouble[1] = Double.NEGATIVE_INFINITY;
		tmpDouble[2] = BAC.forceCpg;
		tmpDouble[3] = 0.0;
		cytosineListMap.put("CG-1".toUpperCase(), tmpDouble);
		tmpDouble = new Double[4];
		tmpDouble[0] = Double.NEGATIVE_INFINITY;
		tmpDouble[1] = Double.NEGATIVE_INFINITY;
		tmpDouble[2] = BAC.forceChg;
		tmpDouble[3] = 0.0;
		cytosineListMap.put("CHG-1".toUpperCase(), tmpDouble);
		tmpDouble = new Double[4];
		tmpDouble[0] = Double.NEGATIVE_INFINITY;
		tmpDouble[1] = Double.NEGATIVE_INFINITY;
		tmpDouble[2] = BAC.forceChh;
		tmpDouble[3] = 0.0;
		cytosineListMap.put("CHH-1".toUpperCase(), tmpDouble);
		if(BAC.sequencingMode == MethylSNPModel.GM){
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceGch;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("GCH-2".toUpperCase(), tmpDouble);
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceGcg;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("GCG-2".toUpperCase(), tmpDouble);
			tmpDouble = new Double[4];
			tmpDouble[0] = Double.NEGATIVE_INFINITY;
			tmpDouble[1] = Double.NEGATIVE_INFINITY;
			tmpDouble[2] = BAC.forceHcg;
			tmpDouble[3] = 0.0;
			cytosineListMap.put("HCG-2".toUpperCase(), tmpDouble);
		}
		if(!BAC.autoEstimateOtherCytosine.isEmpty()){
			String[] tmpArray = BAC.autoEstimateOtherCytosine.split(";");
			for(String tmp : tmpArray){
				tmpDouble = new Double[4];
				tmpDouble[0] = Double.NEGATIVE_INFINITY;
				tmpDouble[1] = Double.NEGATIVE_INFINITY;
				tmpDouble[2] = Double.NaN;
				tmpDouble[3] = 0.0;
				cytosineListMap.put(tmp.toUpperCase(), tmpDouble);
				int tmpLength = tmp.split("-")[0].length();
				if(tmpLength > maxCytosineLength){
					maxCytosineLength = tmpLength;
				}
			}
		}
		
		if(!BAC.forceOtherCytosine.isEmpty()){
			String[] tmpArray = BAC.forceOtherCytosine.split(";");
			for(String tmp : tmpArray){
				String[] mapElement = tmp.split(":");
				tmpDouble = new Double[4];
				tmpDouble[0] = Double.NEGATIVE_INFINITY;
				tmpDouble[1] = Double.NEGATIVE_INFINITY;
				tmpDouble[2] = Double.parseDouble(mapElement[1]);
				tmpDouble[3] = 0.0;
				cytosineListMap.put(mapElement[0].toUpperCase(), tmpDouble);
				int tmpLength = mapElement[0].split("-")[0].length();
				if(tmpLength > maxCytosineLength){
					maxCytosineLength = tmpLength;
				}
			}
		}
	}
	
	public CytosineTypeStatus clone(){
		CytosineTypeStatus cts = new CytosineTypeStatus(BAC);
		cts.isC = this.isC;
		cts.isCpg = this.isCpg;
		cts.isChg = this.isChg;
		cts.isChh = this.isChh;
		cts.isGch = this.isGch;
		cts.isGcg = this.isGcg;
		cts.isHcg = this.isHcg;

		cts.cytosineMethyLevel = this.cytosineMethyLevel;
		cts.cpgMethyLevel = this.cpgMethyLevel;
		cts.chgMethyLevel = this.chgMethyLevel;
		cts.chhMethyLevel = this.chhMethyLevel;
		cts.gchMethyLevel = this.gchMethyLevel;
		cts.gcgMethyLevel = this.gcgMethyLevel;
		cts.hcgMethyLevel = this.hcgMethyLevel;
		
		
		cts.BAC = this.BAC;
		
		cts.maxCytosineLength = this.maxCytosineLength;
		cts.cytosineListMap = new HashMap<String, Double[]>();
		
		for(String key : this.cytosineListMap.keySet()){
			Double[] value = this.cytosineListMap.get(key);
			value[0] = Double.NEGATIVE_INFINITY;
			value[1] = Double.NEGATIVE_INFINITY;
			value[3] = 0.0;
			cts.cytosineListMap.put(key, value);
		}
		
		return cts;
		
	}
	
}
