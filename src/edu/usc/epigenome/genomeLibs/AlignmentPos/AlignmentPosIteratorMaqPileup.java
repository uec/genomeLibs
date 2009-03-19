package edu.usc.epigenome.genomeLibs.AlignmentPos;

import java.io.IOException;

import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.StrandedFeature;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MiscUtils;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosDepthOnly;
import edu.usc.epigenome.genomeLibs.AlignmentPos.AlignmentPosSnps;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPosRich;

public class AlignmentPosIteratorMaqPileup extends AlignmentPosIterator {

	//private int totalBasesRead = 0;
	
	private AlignmentPos hasNextCache = null;
	
	public AlignmentPosIteratorMaqPileup(String fn, AlignmentPosOptions apos) 
	throws IOException {
		super(fn, apos);
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.genomeLibs.AlignmentPosIterator#hasNext()
	 */
	@Override
	public boolean hasNext() {
		
		if (hasNextCache != null) return true;
		
		boolean out = false;
		try
		{
			AlignmentPos nextAp = this.nextAlignment(false);
			hasNextCache = nextAp;
			out = (nextAp != null);
//			AlignmentPos nextAp = this.nextAlignment(true);
//			out =  (nextAp != null);
		}
		catch (Exception e)
		{
			System.err.println("Could not process file " + openFile + "\n" + e.toString());
			e.printStackTrace(System.err);
			System.exit(0);
		}
		return out;
	}

	@Override
	protected AlignmentPos nextAlignment()
	throws IOException, IllegalSymbolException, Exception
	{
		return nextAlignment(false);
	}

	protected AlignmentPos nextAlignment(boolean rollback)
	throws IOException, IllegalSymbolException, Exception
	{
		if (hasNextCache != null)
		{
			AlignmentPos ap = hasNextCache;
			hasNextCache = null;  // Reset the cache 
			return ap;
		}
	
		
		//		String line = this.openStream.readLine();
		//		String[] line_items = line.split("\t");  //TODO VERY SLOW, 25% of execution time

		if (rollback) this.openStream.mark(10000);
		
		AlignmentPos ap = null;
		boolean done = false;

		while (!done && (ap == null))
		{
			String[] line_items = ListUtils.readLineSplitByChar(this.openStream, '\t', 20);
			
			if (line_items == null)
			{
				done = true;
			}
			else if (line_items.length == 0)
			{
				// Blank line, keep trying (unless we hit the end of file)
				done = !this.openStream.ready();
			}
			else if ((line_items.length == 7) || (line_items.length == 8)) // length 7 if no reads
			{
				String line_chr = line_items[0];
				int line_pos = Integer.parseInt(line_items[1]);
				char line_ref = line_items[2].charAt(0);
				//		int line_count = Integer.parseInt(line_items[3]);
				int read_count = Integer.parseInt(line_items[3]);
				String snps_str = line_items[4];
				char[] snps = snps_str.toCharArray();
				String base_quals = line_items[5];
				//		String mapping_quals = line_items[6];
				
				String read_positions = "";
				if ((line_items.length>7) && (read_count>0)) read_positions = line_items[7];


				// Make the output object.  Just make one with SNPs, and then reduce if necessary
				if (this.apOptions.trackBisulfiteConversion)
				{
					ap = new AlignmentPosSnpsBisulfiteConverted(line_ref, line_chr, line_pos, this.apOptions);
				}
				else
				{
					ap = new AlignmentPosSnps(line_ref, line_chr, line_pos, this.apOptions);
				}
				
				ap.setStrand(StrandedFeature.POSITIVE); // Positive by default
				//	System.err.println("ap=" + ap);
				
				if (read_count>0)
				{
				addMaqPositions(this.apOptions, (AlignmentPosSnps)ap, snps, base_quals, read_positions);
				}

				if (!apOptions.trackSnps)
				{
					AlignmentPosDepthOnly newAp = new AlignmentPosDepthOnly(ap);
					newAp.setDepth(ap.getDepth());
					ap = newAp;
					//			System.err.println("ap=" + ap);
				}
			}
			else
			{
				System.err.println("Illegal Maq pileup line: " + ListUtils.excelLine(line_items));
				//throw new Exception("Illegal Maq pileup line: " + ListUtils.excelLine(line_items));
			}
		}
		
		if (rollback) this.openStream.reset();
		return ap;
	}
	
	
	
	protected static void addMaqPositions(AlignmentPosOptions inApOptions, AlignmentPosSnps ap, char[] snps, String baseQualsStr, String readPositionsStr)
	throws IllegalSymbolException
	{
		//String[] readPositionsStrs = readPositionsStr.split(","); //TODO VERY SLOW, 25% of execution time, but it's probably the best we can do with a string
		String[] readPositionsStrs = ListUtils.splitByChar(readPositionsStr, ','); 

		
		char[] baseQualChars = baseQualsStr.toCharArray();
		
		for (int i=0; i < (snps.length-1); i++)
		{
			int qual = MiscUtils.fastqQualCodeToInt(baseQualChars[i+1], false); // First one is a "@" char
			
			if (qual >= inApOptions.minQualityScore)
			{
				char snpChar = snps[i+1]; // First one is a "@" char
				ReadPos rp = maqPileupCharToReadPos(snpChar, ap.ref);

				int cycle = Integer.parseInt(readPositionsStrs[i]);

				if (!inApOptions.onlyFirstCycle || (cycle==1))
				{
					if (inApOptions.trackPositions || inApOptions.trackQuals)
					{
						int apCycle = (inApOptions.trackPositions) ? cycle : ReadPos.UNKNOWN_CYCLE;
						int apQual = (inApOptions.trackQuals) ? qual : ReadPos.UNKNOWN_QUAL;
						rp = new ReadPosRich(rp, apCycle, apQual);
					}
					ap.add(rp);
				}
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

		return new ReadPos(readstrandC,readstrandForwardStrand ? StrandedFeature.POSITIVE : StrandedFeature.NEGATIVE);
	}
	
	

}
