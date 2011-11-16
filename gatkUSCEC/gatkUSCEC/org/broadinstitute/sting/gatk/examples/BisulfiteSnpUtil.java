package org.broadinstitute.sting.gatk.examples;

public class BisulfiteSnpUtil {
	
	public static boolean isCytosine(byte base, boolean bisulfiteConversionSpace)
	{
		char refC = (char) base;
		boolean out;
		
		if (bisulfiteConversionSpace)
		{
			out = ((refC == 'C') || (refC == 'T'));
		}
		else
		{
			out = (refC == 'C');
		}
		
		return out; 
	}
	
	public static boolean isCytosine(int pos, String seqStr, boolean bisulfiteConversionSpace)
	{
		char refC = seqStr.charAt(pos);
		
		boolean out;
		
		if (bisulfiteConversionSpace)
		{
			out = ((refC == 'C') || (refC == 'T'));
		}
		else
		{
			out = (refC == 'C');
		}
		
		return out; 
	}
	
	public static boolean isThymine(int pos, String refStr)
	{
		char refC = refStr.charAt(pos);
		
		return (refC == 'T'); 
	}
	
}
