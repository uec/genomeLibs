package edu.usc.epigenome.genomeLibs;

import org.broadinstitute.sting.utils.BaseUtils;

public class GatkBaseUtils {
	public static char CharFromBaseByte(byte inBase)
	{
		return BaseUtils.baseIndexToSimpleBaseAsChar(BaseUtils.extendedBaseToBaseIndex(inBase));
	}
	
	public static byte BaseByteFromChar(char inBaseChar)
	{
		byte out = BaseUtils.N;
		switch (inBaseChar)
		{
		case 'A':
		case 'a':
			out = BaseUtils.A;
			break;
		case 'C':
		case 'c':
			out = BaseUtils.C;
			break;
		case 'T':
		case 't':
			out = BaseUtils.T;
			break;
		case 'G':
		case 'g':
			out = BaseUtils.G;
			break;
		case 'N':
		case 'n':
		case '0':
		default:
			out = BaseUtils.N;
			break;
		}
		
		return out;
	}
}
