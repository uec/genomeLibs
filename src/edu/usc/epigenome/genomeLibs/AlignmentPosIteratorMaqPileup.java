package edu.usc.epigenome.genomeLibs;

import java.io.IOException;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

public class AlignmentPosIteratorMaqPileup extends AlignmentPosIterator {

	//private int totalBasesRead = 0;
	
	public AlignmentPosIteratorMaqPileup(String fn, AlignmentPosOptions apos) 
	throws IOException {
		super(fn, apos);
	}

	@Override
	protected AlignmentPos nextAlignment()
	throws IOException, IllegalSymbolException
	{
//		String line = this.openStream.readLine();
//		String[] line_items = line.split("\t");  //TODO VERY SLOW, 25% of execution time

		String[] line_items = ListUtils.readLineSplitByChar(this.openStream, '\t', 20);
		
		String line_chr = line_items[0];
		int line_pos = Integer.parseInt(line_items[1]);
		char line_ref = line_items[2].charAt(0);
//		int line_count = Integer.parseInt(line_items[3]);
		String snps_str = line_items[4];
		char[] snps = snps_str.toCharArray();
		String base_quals = line_items[5];
//		String mapping_quals = line_items[6];
		String read_positions = line_items[7];
		
		
		// Make the output object.  Just make one with SNPs, and then reduce if necessary
		AlignmentPos ap = new AlignmentPosSnps(line_ref, line_chr, line_pos, this.apOptions);
//		System.err.println("ap=" + ap);
		this.addMaqPositions((AlignmentPosSnps)ap, snps, base_quals, read_positions);

		if (!apOptions.trackSnps)
		{
			AlignmentPosDepthOnly newAp = new AlignmentPosDepthOnly(ap);
			newAp.setDepth(ap.getDepth());
			ap = newAp;
//			System.err.println("ap=" + ap);
		}
		
		return ap;
	}
	
	
	
	protected void addMaqPositions(AlignmentPosSnps ap, char[] snps, String baseQualsStr, String readPositionsStr)
	throws IllegalSymbolException
	{
		//String[] readPositionsStrs = readPositionsStr.split(","); //TODO VERY SLOW, 25% of execution time, but it's probably the best we can do with a string
		String[] readPositionsStrs = ListUtils.splitByChar(readPositionsStr, ','); 

		
		char[] baseQualChars = baseQualsStr.toCharArray();
		
		for (int i=0; i < (snps.length-1); i++)
		{
			int qual = fastqQualCodeToInt(baseQualChars[i+1]); // First one is a "@" char
			
			if (qual >= apOptions.minQualityScore)
			{
				char snpChar = snps[i+1]; // First one is a "@" char
				ReadPos rp = maqPileupCharToReadPos(snpChar, ap.ref);

				if (apOptions.trackPositionsQuals)
				{
					int readPos = Integer.parseInt(readPositionsStrs[i]);
					rp = new ReadPosRich(rp, readPos, qual);
				}
				
				ap.add(rp);
			}
			
		}
		
	}
	
	protected static ReadPos maqPileupCharToReadPos(char c, Symbol ref)
	throws IllegalSymbolException
	{
		// Check for maq special characters
		Symbol readstrandC;
		Boolean readstrandForwardStrand;
		
		switch(c)
		{
		case ',': readstrandC = ref; readstrandForwardStrand = true; break;
		case '.': readstrandC = ref; readstrandForwardStrand = false; break;
		default: readstrandC = DNATools.forSymbol(c); readstrandForwardStrand = Character.isUpperCase(c); break;
		}

		return new ReadPos(readstrandC,readstrandForwardStrand);
	}
	
	
	public static int fastqQualCodeToInt(char c)
	{
		//int v = Character.getNumericValue(c) - 33;
		int v = (int)c - 33;
		return v;
	}
	

}
