package edu.usc.epigenome.genomeLibs;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

/**
 * @author benb
 *
 */
/**
 * @author benb
 *
 */
public class PicardUtils {

	static final int BUFFERLEN = 1000;
	
	
	static Pattern mdPat = Pattern.compile("\\G(?:([0-9]+)|([ACTGNactgn])|(\\^[ACTGNactgn]+))");
	
	/*
	 * ROLL BACK INTO PICARD
	 */
	

	public static String getReadString(SAMRecord rec, boolean relativeToReadStrand)
	{
		String seq = rec.getReadString().toUpperCase();
		if (relativeToReadStrand && rec.getReadNegativeStrandFlag())
		{
			seq = MiscUtils.revCompNucStr(seq);
		}
		return seq;
	}

	public static String getBaseQualityString(SAMRecord rec, boolean relativeToReadStrand)
	{
		String baseQual = rec.getBaseQualityString();
		if (relativeToReadStrand && rec.getReadNegativeStrandFlag())
		{
			baseQual = MiscUtils.revString(baseQual);
		}
		return baseQual;
	}

	public static String refStr(SAMRecord rec, boolean relativeToReadStrand)
	throws Exception
	{
		String md = (String)rec.getAttribute("MD");
		Cigar cigar = rec.getCigar();
		String seq = rec.getReadString().toUpperCase();
		String ref = PicardUtils.refStrFromMd(seq, md, cigar).toUpperCase();
		if (relativeToReadStrand && rec.getReadNegativeStrandFlag())
		{
			ref = MiscUtils.revCompNucStr(ref);
		}
		return ref;
	}

	public static String refStrFromMd(String seq, String md, Cigar cigar)
	throws Exception
	{
		if (seq == null) throw new Exception("Can not run refStrFromMd with a null seq variable");
		if (md == null) throw new Exception("Can not run refStrFromMd with a null seq variable");
		
		// Use sb as the reference output string
		StringBuilder sb = new StringBuilder(BUFFERLEN);

		Matcher match = mdPat.matcher(md);
		int curSeqPos = 0;
		//int curMdPos = 0; // Not the same as seq pos when you have indels

		int savedBases = 0;
		for (final CigarElement cigEl : cigar.getCigarElements()) 
		{
			int cigElLen = cigEl.getLength();
			CigarOperator cigElOp = cigEl.getOperator();
//			System.err.printf("\tCigar El: len=%d, op=%s, consumesRead=%b, consumesRef=%b\n",
//					cigElLen,cigElOp,cigElOp.consumesReadBases(), cigElOp.consumesReferenceBases());
			
			
			// If it consumes reference bases, it's either a match or a deletion in the sequence
			// read.  Either way, we're going to need to parse throught the MD.
			if (cigElOp.consumesReferenceBases())
			{
				// We have a match region, go through the MD
				int basesMatched = 0;
				
				// Do we have any saved matched bases?
				while ((savedBases>0) && (basesMatched < cigElLen))
				{
					sb.append(seq.charAt(curSeqPos++));
					savedBases--;
					basesMatched++;
//					System.err.printf("\t\tDepleting saved bases, saved=%d, curSeqPos=%d, basesMatched=%d\n",savedBases,curSeqPos,basesMatched); 
				}

				while (basesMatched < cigElLen)
				{
					boolean matched = match.find();
					if (matched)
					{
//						System.err.println("Matched , basesMatched=" + basesMatched + ", match=" + match.group() + "," + match.group(1) + "," + match.group(2) + ", start=" + match.start());
						String mg;
						if ( ((mg = match.group(1)) !=null) && (mg.length() > 0) )
						{
							// It's a number , meaning a series of matches
							int num = Integer.parseInt(mg);
							for (int i = 0; i < num; i++)
							{
								if (basesMatched<cigElLen)
								{
									sb.append(seq.charAt(curSeqPos));
									curSeqPos++;
								}
								else
								{
									savedBases++;
								}
								basesMatched++;
							}
						}

						else if ( ((mg = match.group(2)) !=null) && (mg.length() > 0) )
						{
							// It's a single nucleotide, meaning a mismatch
							if (basesMatched<cigElLen)
							{
								sb.append(mg.charAt(0));
								curSeqPos++;
							}
							else
							{
								savedBases++;
							}
							basesMatched++;
						}
						else if ( ((mg = match.group(3)) !=null) && (mg.length() > 0) )
						{
							// It's a deletion, starting with a caret
							// don't include caret
							for (int i = 1; i < mg.length(); i++)
							{
								// Since this function is actually just meant to make a reference that lines up nucleotide 
								//  for nucleotide with the sequence read, we don't actually add the insertion to the reference.
								//sb.append(mg.charAt(i));
								basesMatched++;
							}
							
							// Check just to make sure.
							if (basesMatched != cigElLen)
							{
								throw new Exception("Got a deletion in CIGAR (" + cigar + ", deletion " + cigElLen + 
										" length) with an unequal ref insertion in MD (" + md + ", md " + basesMatched + " length");
							}
							if (cigElOp != CigarOperator.DELETION)
							{
								throw new Exception ("Got an insertion in MD ("+md+") without a corresponding deletion in cigar ("+cigar+")");
							}
							
						}
						else
						{
							matched = false;
						}
					}

					if (!matched)
					{
						throw new Exception("Illegal MD pattern: " + md);
					}

//					System.err.println("SavedBasesMatched=" + savedBases);
				}

			}
			else if (cigElOp.consumesReadBases())
			{
				// We have an insertion in read
				for (int i = 0; i < cigElLen; i++)
				{
					char c = (cigElOp == CigarOperator.SOFT_CLIP) ? '0' : '-';
					sb.append( c );
					curSeqPos++;
				}
			}
			else
			{
				// It's an op that consumes neither read nor reference bases.  Do we just ignore??
			}

		}
		
		return sb.toString();
	}

	/**
	 * @param samRecord
	 * @return true if the read or the reference sequence contains at least one CpG
	 */
	public static boolean readContainsCpg(SAMRecord samRecord) throws Exception {
		int nCpgs = readNumCpgs(samRecord);
		return (nCpgs>=1);
	}

	/**
	 * @param samRecord
	 * @return Returns the number of CpGs in the sequence, either in the read or reference.
	 */
	private static int readNumCpgs(SAMRecord samRecord) throws Exception {
		String ref = PicardUtils.refStr(samRecord, true);

		int out = 0;
		for (int i = 0; i < ref.length(); i++)
		{
			if (isCpg(i,ref)) out++;
		}
		//System.err.printf("\t%d cpgs\n",out);
		
		return out;
	}
		
	public static char nextBaseSeq(SAMRecord samRecord, int pos)
	{
		return nextBaseSeq(samRecord, pos, 20);
	}

	
	/**
	 * @param samRecord
	 * @param pos IMPORTANT: must be relative to the 5' end of the read, not the genome assembly
	 * @param minBaseQual
	 * @return
	 */
	public static char nextBaseSeq(SAMRecord samRecord, int pos, int minBaseQual)
	{
		String seq = PicardUtils.getReadString(samRecord, true);
		byte[] quals = samRecord.getBaseQualities();
		return nextBaseSeq(pos, seq, false, quals, minBaseQual);
	}
	


	public static char nextBaseSeq(int pos, String seqStr, byte[] quals, int minBaseQual)
	{
		return nextBaseSeq(pos, seqStr, false, quals, minBaseQual);
	}

	public static char nextBaseSeq(int pos, String seqStr)
	{
		return nextBaseSeq(pos, seqStr, false);
	}

	public static char nextBaseSeq(int pos, String seqStr, boolean revStrand)
	{
		return nextBaseSeq(pos, seqStr, revStrand, null, 0);
	}

	public static char nextBaseSeq(int pos, String seqStr, boolean revStrand, byte[] quals, int minBaseQual)
	{
		int prePos = (revStrand) ? (pos-1) : (pos+1);
		char out = '0';
		if ((prePos>=0) && (prePos<seqStr.length()))
		{
			if ((quals != null) && (quals[prePos]<minBaseQual))
			{
				out = 'N';
			}
			else
			{
				out = seqStr.charAt(prePos);
				if (revStrand) out = MiscUtils.revCompNuc(out);
			}
		}
		return out;		
	}
	
	public static char nextBaseRef(int pos, String refStr)
	{
		return nextBaseRef(pos, refStr, false);
	}
	
	/**
	 * @param pos  This is relative to refStr strand regardless of value of revStrand
	 * @param refStr
	 * @param revStrand  If this is true, we give the prior base, reverse complemented.
	 * @return
	 */
	public static char nextBaseRef(int pos, String refStr, boolean revStrand)
	{
		char refCnext = '0';

		if (revStrand) 
		{
			if (pos == 0) return '0'; // At the last character
			refCnext = MiscUtils.revCompNuc(refStr.charAt(pos-1));
		}
		else
		{
			if (pos >= (refStr.length()-1)) return '0'; // At the last character
			refCnext = refStr.charAt(pos+1);
		}
		
		return refCnext;
	}
	
	public static char preBaseSeq(SAMRecord samRecord, int pos)
	{
		return preBaseSeq(samRecord, pos, 20);
	}
	
	/**
	 * @param samRecord
	 * @param pos IMPORTANT: Must be relative to the 5' end of the read, not the genome assembly
	 * @param minBaseQual
	 * @return
	 */
	public static char preBaseSeq(SAMRecord samRecord, int pos, int minBaseQual)
	{
		String seq = PicardUtils.getReadString(samRecord, true);
		byte[] quals = samRecord.getBaseQualities();
		return preBaseSeq(pos, seq, false, quals, minBaseQual);
	}
	
	public static char preBaseSeq(int pos, String seqStr)
	{
		return preBaseSeq(pos, seqStr, false);
	}

	public static char preBaseSeq(int pos, String seqStr, boolean revStrand)
	{
		return preBaseSeq(pos, seqStr, revStrand, null, 0);
	}
	
	public static char preBaseSeq(int pos, String seqStr, byte[] quals, int minBaseQual)
	{
		return preBaseSeq(pos, seqStr, false, quals, minBaseQual);
	}
	
	public static char preBaseSeq(int pos, String seqStr, boolean revStrand, byte[] quals, int minBaseQual)
	{
		int prePos = (revStrand) ? (pos+1) : (pos-1);
		char out = '0';
		if ((prePos>=0) && (prePos<seqStr.length()))
		{
			if ((quals != null) && quals[prePos]<minBaseQual)
			{
				out = 'N';
			}
			else
			{
				out = seqStr.charAt(prePos);
				if (revStrand) out = MiscUtils.revCompNuc(out);
//				System.err.printf("\tprebase seq (rev=%s) out = %c\n",revStrand,out);
			}
		}
		return out;		
	}
	
	public static char preBaseRef(int pos, String refStr)
	{
		return preBaseRef(pos, refStr, false);
	}
	
	/**
	 * @param pos  This is relative to refStr strand regardless of value of revStrand
	 * @param refStr
	 * @param revStrand  If this is true, we give the prior base, reverse complemented.
	 * @return
	 */
	public static char preBaseRef(int pos, String refStr, boolean revStrand)
	{
		char refCpre = '0';

		if (revStrand) 
		{
			if (pos >= refStr.length()-1) return '0'; // At the last character
			refCpre = MiscUtils.revCompNuc(refStr.charAt(pos+1));
		}
		else
		{
			if (pos == 0) return '0'; // At the last character
			refCpre = refStr.charAt(pos-1);
		}
		
		return refCpre;
	}
	
	
	
	public static boolean isCpg(int pos, String refStr)
	{
		if (pos >= (refStr.length()-1)) return false; // At the last character
		
		char refCnext = refStr.charAt(pos+1);
		
		return ( isCytosine(pos,refStr,false) && (refCnext == 'G') );
	}	

	
	/**
	 * @param pos
	 * @param refStr
	 * @param seqStr If provided, we count if either the reference or the seq is a CpG
	 * @return
	 */
	public static boolean isCpg(int pos, String refStr, String seqStr)
	{
		if (pos >= (refStr.length()-1)) return false; // At the last character
		
		char refCnext = refStr.charAt(pos+1);
		char seqCnext = seqStr.charAt(pos+1);
		boolean out = isCytosine(pos,refStr,false) && ((refCnext == 'G') || (seqCnext == 'G'));
		return out;
	}	

	// The G opposit the CpG
	public static boolean isOppositeCpg(int pos, String refStr)
	{
		if (pos == 0) return false; // At the first char
		
		char refCprev = refStr.charAt(pos-1);
		
		return ( isGuanine(pos,refStr) && (refCprev == 'C') );
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
	
	public static boolean isGuanine(int pos, String refStr)
	{
		char refC = refStr.charAt(pos);
		
		return (refC == 'G') ; 
	}
	
	public static boolean isAdenine(int pos, String refStr)
	{
		char refC = refStr.charAt(pos);
		
		return (refC == 'A') ; 
	}
	
	public static char revNucleotide(char nucleotide)
	throws Exception
	{
		switch(nucleotide){
		case 'N':
			break;
		case 'A':
			nucleotide = 'T';
			break;
		case 'T':
			nucleotide = 'A';
			break;
		case 'G':
			nucleotide = 'C';
			break;
		case 'C':
			nucleotide = 'G';
			break;
		default:
			throw new Exception("Can't recognize seq char: " + nucleotide);
		}
		return nucleotide;
		
	}
	
	
	public static boolean isConverted(int pos, String refStr, String seqStr)
	{
		char refC = refStr.charAt(pos);
		char seqC = seqStr.charAt(pos);
		
		return ((refC == 'C') && (seqC == 'T'));
	}

	public static boolean isGch(int pos, String refStr)
	{
		if (pos <= 0) return false;
		if (pos == (refStr.length()-1)) {
			if (isGuanine(pos-1,refStr) && isCytosine(pos,refStr,false)) {
				return true;
			}
			else{
				return false;
			}
		} // At the last character
		if (pos > (refStr.length()-1)) return false;
		return ( isGuanine(pos-1,refStr) && isCytosine(pos,refStr,false) && !isGuanine(pos+1,refStr) );
	}	
	
	public static boolean isOppositeGch(int pos, String refStr)
	{
		if (pos >= refStr.length()-1) return false;
		if (pos == 0) {
			if (isGuanine(pos,refStr) && isCytosine(pos+1,refStr,false)) {
				return true;
			}
			else{
				return false;
			}
		} // At the first character
		if (pos < 0) return false;
		return ( isCytosine(pos+1,refStr,false) && isGuanine(pos,refStr) && !isCytosine(pos-1,refStr,false) );
	}	

/*	public static boolean isOppositeGch(int pos, String refStr){
		if (pos > (refStr.length()-2)) return false;
		if (pos == (refStr.length()-2)) return (isCytosine(pos+1,refStr,false) && isGuanine(pos,refStr));
		return (isCytosine(pos+1,refStr,false) && isGuanine(pos,refStr) && !isGuanine(pos+2,refStr));
	}
*/	
	
	public static boolean isGcg(int pos, String refStr)
	{
		if (pos >= (refStr.length()-1)) return false; // At the last character
		if (pos <= 0) return false;
		
		return ( isGuanine(pos-1,refStr) && isCytosine(pos,refStr,false) && isGuanine(pos+1,refStr) );
	}	
	
	public static boolean isOppositeGcg(int pos, String refStr)
	{
		if (pos >= (refStr.length()-1)) return false; // At the last character
		if (pos <= 0) return false;
		
		return ( isCytosine(pos-1,refStr,false) && isGuanine(pos,refStr) && isCytosine(pos+1,refStr,false) );
	}

/*	public static boolean isOppositeGcg(int pos, String refStr){
		if(pos >= (refStr.length()-2)) return false;
		return isGuanine(pos,refStr) && isCytosine(pos+1,refStr,false) && isGuanine(pos+2,refStr);
	}
*/
	public static boolean isHcg(int pos, String refStr)
	{
		if (pos >= (refStr.length()-1)) return false; // At the last character
		if (pos == 0) {
			if (isGuanine(pos+1,refStr) && isCytosine(pos,refStr,false)) {
				return true;
			}
			else{
				return false;
			}
		}
		
		return ( !isGuanine(pos-1,refStr) && isCytosine(pos,refStr,false) && isGuanine(pos+1,refStr) );
	}	

	
	public static boolean isOppositeHcg(int pos, String refStr)
	{
		if (pos >= (refStr.length()-1)) return false; // At the last character
		if (pos == 0) {
			if (isGuanine(pos,refStr) && isCytosine(pos+1,refStr,false)) {
				return true;
			}
			else{
				return false;
			}
		}
		
		return ( !isCytosine(pos+1,refStr,false) && isGuanine(pos,refStr) && isCytosine(pos-1,refStr,false) );
	}

/*	public static boolean isOppositeHcg(int pos, String refStr){
		if(pos <= 0) return false;
		if(pos == 1) return (isGuanine(pos,refStr) && isCytosine(pos-1,refStr,false));
		return ( !isGuanine(pos-2,refStr) && isGuanine(pos,refStr) && isCytosine(pos-1,refStr,false) );
		
	}
*/	
	
	public static boolean isHch(int pos, String refStr)
	{
		if (pos > (refStr.length()-1) || pos < 0) return false; // At the last character
		if (pos == 0){
			if (isCytosine(pos,refStr,false) && !isGuanine(pos+1,refStr)){
				return true;
			}
			else{
				return false;
			}
		}
		if (pos == refStr.length()-1){
			if (!isGuanine(pos-1,refStr) && isCytosine(pos,refStr,false)){
				return true;
			}
			else{
				return false;
			}
		}
		
		return ( !isGuanine(pos-1,refStr) && isCytosine(pos,refStr,false) && !isGuanine(pos+1,refStr) );
	}	
	
	public static boolean isOppositeHch(int pos, String refStr)
	{
		if (pos > (refStr.length()-1) || pos < 0) return false;
		if (pos == 0){
			if (!isCytosine(pos+1,refStr,false) && isGuanine(pos,refStr)){
				return true;
			}
			else{
				return false;
			}
		}
		if (pos == refStr.length()-1){
			if (isGuanine(pos,refStr) && !isCytosine(pos-1,refStr,false)){
				return true;
			}
			else{
				return false;
			}
		}
		
		return ( !isCytosine(pos-1,refStr,false) && isGuanine(pos,refStr) && !isCytosine(pos+1,refStr,false) );
	}

	
}
