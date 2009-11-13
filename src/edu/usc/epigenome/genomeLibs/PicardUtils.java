package edu.usc.epigenome.genomeLibs;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;

public class PicardUtils {

	static final int BUFFERLEN = 1000;
	
	
	static Pattern pairPat = Pattern.compile("\\G(?:([0-9]+)|([ACTGNactgn])|(\\^[ACTGNactgn]+))");
	
	/*
	 * ROLL BACK INTO PICARD
	 */

	public static String refStrFromMd(String seq, String md, Cigar cigar)
	throws Exception
	{
		// Use sb as the reference output string
		StringBuilder sb = new StringBuilder(BUFFERLEN);

		Matcher match = pairPat.matcher(md);
		int curSeqPos = 0;
		//int curMdPos = 0; // Not the same as seq pos when you have indels

		int savedBases = 0;
		for (final CigarElement cigEl : cigar.getCigarElements()) 
		{
			int cigElLen = cigEl.getLength();
			CigarOperator cigElOp = cigEl.getOperator();
			System.err.printf("\tCigar El: len=%d, op=%s, consumesRead=%b, consumesRef=%b\n",
					cigElLen,cigElOp,cigElOp.consumesReadBases(), cigElOp.consumesReferenceBases());
			
			
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
					System.err.printf("\t\tDepleting saved bases, saved=%d, curSeqPos=%d, basesMatched=%d\n",savedBases,curSeqPos,basesMatched); 
				}

				while (basesMatched < cigElLen)
				{
					boolean matched = match.find();
					if (matched)
					{
						System.err.println("Matched , basesMatched=" + basesMatched + ", match=" + match.group() + "," + match.group(1) + "," + match.group(2) + ", start=" + match.start());
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

					System.err.println("SavedBasesMatched=" + savedBases);
				}

			}
			else if (cigElOp.consumesReadBases())
			{
				// We have an insertion in read
				for (int i = 0; i < cigElLen; i++)
				{
					sb.append( (cigElOp.consumesReadBases() && !cigElOp.consumesReferenceBases()) ? "-" : seq.charAt(curSeqPos) );
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
		
	
	
	
}
