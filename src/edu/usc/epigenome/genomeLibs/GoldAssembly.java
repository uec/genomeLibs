package edu.usc.epigenome.genomeLibs;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

import org.biojava.bio.seq.Sequence;
import org.biojavax.bio.seq.RichSequence.IOTools;
import org.biojavax.bio.seq.RichSequenceIterator;


public class GoldAssembly {

	public static Map<String,Integer> c_chr_map; // "mm9_chr1"

	
	public static int chromLengthStatic(String chr, String genome)
	throws Exception
	{
		int len = 0;
		
		String key = genome + "__" + chr;
		Integer len_int = c_chr_map.get(key);
		if	(len_int == null)
		{
			System.err.println("Number of keys in c_chr_map = " + c_chr_map.size());
			throw new Exception("Can not identify chrom lenth for " + key);
		}
		else
		{
			len = len_int.intValue();
		}

		System.err.println("Finding length of chrom " + chr + 
				" from " + key + ", len=" + len);
		
		return len;
	}
	
	
	
	/***
	 * 
	 * @param genome
	 * @param extended: Extended contains all chromosomes containing underscore .. 
	 * randoms, haplotypes, etc.  Not recommended for general re-sequencing analysis 
	 * @return
	 */
	public static Iterator<String> chromIterator(String genome, boolean extended)
	{
		Iterator<String> map_it = c_chr_map.keySet().iterator();
		
		Set<String> out_set = new TreeSet<String>();
		while (map_it.hasNext())
		{
			String full = map_it.next();
			String[] parts = full.split("__");
			
			if (parts[0].equals(genome))
			{
				String rest = parts[1]; // Only valid if chrom can not contain __ delimiter
				
				if (extended || (!rest.contains("_")))
				{
					out_set.add(rest);
				}
			}
		}
		
		return out_set.iterator();
	}
		
	
	/************** FETCHING SEQUENCES
	 * 
	 */

	public static Sequence chromSeq(String genome, String chr)
	throws Exception
	{

		String suffix = "fa";
		if (genome.matches("^tair.*"))
		{
			suffix = "fas";
		}

		String fn = genomeDataDir(genome) + "/chromosomes/" + chr + "." + suffix;
		System.err.println("Fetching seq: " + chr);

		BufferedReader br;
		try
		{
			br = new BufferedReader(new FileReader(fn));
		}
		catch (Exception e)
		{
			System.err.println("Could not find genomic fasta file: " + fn);
			return null;
		}

		RichSequenceIterator stream = IOTools.readFastaDNA(br, null);
		Sequence seq = stream.nextRichSequence();

		return seq;
	}

	public static String genomeDataDir(String genome)
	{

		if (genome.endsWith("tair7"))
		{
			// They are the same, except for a small number of SNPs
			genome = "tair8";
		}

		return "/Users/benb/genome_data/ucsc/goldenPath/" + genome;
	}


	// Static init
	static {
		c_chr_map = new HashMap<String,Integer>();

		// Got from UCSC
		//
		// Transformed with:
		// perl -n -e '@a=split(/\s+/); $chr=$a[0]; $l=$a[1]; $l=~s/,//g; print "c_chr_map.put(\"$chr\", new Integer($l));\n";'

		// Human
		c_chr_map.put("hg18__chr1", new Integer(247249719));
		c_chr_map.put("hg18__chr1_random", new Integer(1663265));
		c_chr_map.put("hg18__chr2", new Integer(242951149));
		c_chr_map.put("hg18__chr2_random", new Integer(185571));
		c_chr_map.put("hg18__chr3", new Integer(199501827));
		c_chr_map.put("hg18__chr3_random", new Integer(749256));
		c_chr_map.put("hg18__chr4", new Integer(191273063));
		c_chr_map.put("hg18__chr4_random", new Integer(842648));
		c_chr_map.put("hg18__chr5", new Integer(180857866));
		c_chr_map.put("hg18__chr5_h2_hap1", new Integer(1794870));
		c_chr_map.put("hg18__chr5_random", new Integer(143687));
		c_chr_map.put("hg18__chr6", new Integer(170899992));
		c_chr_map.put("hg18__chr6_cox_hap1", new Integer(4731698));
		c_chr_map.put("hg18__chr6_qbl_hap2", new Integer(4565931));
		c_chr_map.put("hg18__chr6_random", new Integer(1875562));
		c_chr_map.put("hg18__chr7", new Integer(158821424));
		c_chr_map.put("hg18__chr7_random", new Integer(549659));
		c_chr_map.put("hg18__chr8", new Integer(146274826));
		c_chr_map.put("hg18__chr8_random", new Integer(943810));
		c_chr_map.put("hg18__chr9", new Integer(140273252));
		c_chr_map.put("hg18__chr9_random", new Integer(1146434));
		c_chr_map.put("hg18__chr10", new Integer(135374737));
		c_chr_map.put("hg18__chr10_random", new Integer(113275));
		c_chr_map.put("hg18__chr11", new Integer(134452384));
		c_chr_map.put("hg18__chr11_random", new Integer(215294));
		c_chr_map.put("hg18__chr12", new Integer(132349534));
		c_chr_map.put("hg18__chr13", new Integer(114142980));
		c_chr_map.put("hg18__chr13_random", new Integer(186858));
		c_chr_map.put("hg18__chr14", new Integer(106368585));
		c_chr_map.put("hg18__chr15", new Integer(100338915));
		c_chr_map.put("hg18__chr15_random", new Integer(784346));
		c_chr_map.put("hg18__chr16", new Integer(88827254));
		c_chr_map.put("hg18__chr16_random", new Integer(105485));
		c_chr_map.put("hg18__chr17", new Integer(78774742));
		c_chr_map.put("hg18__chr17_random", new Integer(2617613));
		c_chr_map.put("hg18__chr18", new Integer(76117153));
		c_chr_map.put("hg18__chr18_random", new Integer(4262));
		c_chr_map.put("hg18__chr19", new Integer(63811651));
		c_chr_map.put("hg18__chr19_random", new Integer(301858));
		c_chr_map.put("hg18__chr20", new Integer(62435964));
		c_chr_map.put("hg18__chr21", new Integer(46944323));
		c_chr_map.put("hg18__chr21_random", new Integer(1679693));
		c_chr_map.put("hg18__chr22", new Integer(49691432));
		c_chr_map.put("hg18__chr22_h2_hap1", new Integer(63661));
		c_chr_map.put("hg18__chr22_random", new Integer(257318));
		c_chr_map.put("hg18__chrx", new Integer(154913754));
		c_chr_map.put("hg18__chrx_random", new Integer(1719168));
		c_chr_map.put("hg18__chry", new Integer(57772954));
		c_chr_map.put("hg18__chrm", new Integer(16571));
		
		// mm8
		c_chr_map.put("mm8__chr1", new Integer(197069962));
		c_chr_map.put("mm8__chr1_random", new Integer(172274));
		c_chr_map.put("mm8__chr2", new Integer(181976762));
		c_chr_map.put("mm8__chr3", new Integer(159872112));
		c_chr_map.put("mm8__chr4", new Integer(155029701));
		c_chr_map.put("mm8__chr5", new Integer(152003063));
		c_chr_map.put("mm8__chr5_random", new Integer(2921247));
		c_chr_map.put("mm8__chr6", new Integer(149525685));
		c_chr_map.put("mm8__chr7", new Integer(145134094));
		c_chr_map.put("mm8__chr7_random", new Integer(243910));
		c_chr_map.put("mm8__chr8", new Integer(132085098));
		c_chr_map.put("mm8__chr8_random", new Integer(206961));
		c_chr_map.put("mm8__chr9", new Integer(124000669));
		c_chr_map.put("mm8__chr9_random", new Integer(17232));
		c_chr_map.put("mm8__chr10", new Integer(129959148));
		c_chr_map.put("mm8__chr10_random", new Integer(10781));
		c_chr_map.put("mm8__chr11", new Integer(121798632));
		c_chr_map.put("mm8__chr12", new Integer(120463159));
		c_chr_map.put("mm8__chr13", new Integer(120614378));
		c_chr_map.put("mm8__chr13_random", new Integer(436191));
		c_chr_map.put("mm8__chr14", new Integer(123978870));
		c_chr_map.put("mm8__chr15", new Integer(103492577));
		c_chr_map.put("mm8__chr15_random", new Integer(105932));
		c_chr_map.put("mm8__chr16", new Integer(98252459));
		c_chr_map.put("mm8__chr17", new Integer(95177420));
		c_chr_map.put("mm8__chr17_random", new Integer(89091));
		c_chr_map.put("mm8__chr18", new Integer(90736837));
		c_chr_map.put("mm8__chr19", new Integer(61321190));
		c_chr_map.put("mm8__chrx", new Integer(165556469));
		c_chr_map.put("mm8__chrx_random", new Integer(39696));
		c_chr_map.put("mm8__chry", new Integer(16029404));
		c_chr_map.put("mm8__chry_random", new Integer(14577732));
		c_chr_map.put("mm8__chrun_random", new Integer(1540053));
		c_chr_map.put("mm8__chrm", new Integer(16299));
	
		// TAIR7
		c_chr_map.put("tair7__chr1", new Integer(30432563));
		c_chr_map.put("tair7__chr2", new Integer(19705359));
		c_chr_map.put("tair7__chr3", new Integer(23470805));
		c_chr_map.put("tair7__chr4", new Integer(18585042));
		c_chr_map.put("tair7__chr5", new Integer(26992728));
		c_chr_map.put("tair7__chrc", new Integer(154478));
		c_chr_map.put("tair7__chrm", new Integer(366924));
	}
	


		


}
