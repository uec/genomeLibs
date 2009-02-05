package edu.usc.epigenome.genomeLibs.ChromScores;

import org.biojava.bio.seq.Sequence;

import edu.usc.epigenome.genomeLibs.GoldAssembly;

public class ChromScoresCpGDensity extends ChromScoresArray {

	public ChromScoresCpGDensity(String genome) 
	throws Exception
	{
		super(genome);
		System.err.println("Initing CpGDensities with genome=" + genome);
	}

	
	public void populateFromFile(String chr)
	throws Exception
	{
		String genome = this.f_genome;
		String fn = "/Users/benb/research/genomic-data-misc/CpG_density/" + 
		chr + "_" + genome + "_cpgdens_200bp.wig";
		System.err.println("Populating densities from " + fn);
		
		this.populateFromFile(chr, fn);
	}

	public void populateFromFile(String chr, String fn)
	throws Exception
	{
		this.populateFromWig(fn);
	}

	
	
	public void populate(String target_chr, int wind)
	throws Exception
	{
		populate(target_chr, wind, 0, 0);
	}

	
	
	// If end_coord == 0 , we do the whole chromosome 
	public void populate(String target_chr, int wind, int start_coord, int end_coord)
	throws Exception
	{
		// Get the chromosome
		Sequence seq = GoldAssembly.chromSeq(f_genome, target_chr);
		String full_residues = seq.seqString().toUpperCase();
		if (end_coord > 0)
		{
			full_residues = full_residues.substring(start_coord, end_coord);
		}	

		populate(target_chr, wind, full_residues, start_coord);
	}
	
	public void populate(String target_chr, int wind, String full_residues, int start_coord)
	throws Exception
	{
		int n = wind+2; // Because our method actually shrinks the window by two
		int wind_half = Math.round((float)wind/(float)2.0);

		
		// Only need to do one strand for CpGs (symmetric)
		int len = full_residues.length();
		int end = len-n;
		int wind_c = 0;
		int wind_g = 0;
		int wind_cg = 0;
		// Must start negative to pick up first window
		for (int i = -n+1; i <= end; i++)
		{
			
			int to = i+n-1;

			char first_one = 'N', second_one = 'N';
			if (i >= 0)
			{
				first_one = full_residues.charAt(i);
				second_one = full_residues.charAt(i+1);
			}
			
			char penultimate_one = 'N';
			if (to >= 1) penultimate_one = full_residues.charAt(to-1); // Unless we are on the very first one
			char last_one = full_residues.charAt(to);
			
			// Add on
			if (last_one == 'G')
			{
				wind_g++;
				if (penultimate_one == 'C') wind_cg++;
			}
			else if (last_one == 'C')
			{
				wind_c++;
			}
			
			// Take off
			if (first_one == 'C')
			{
				wind_c--;
				if (second_one == 'G') wind_cg--;
			}
			else if (first_one == 'G')
			{
				wind_g--;
			}
			
			double c_frac = (double)wind_c/(double)wind;
			double g_frac = (double)wind_g/(double)wind;
			double exp_cg = (double)wind*c_frac*g_frac;
			double density = (exp_cg==0.0) ? 0.0 : (double)wind_cg/exp_cg;
							
			int middle = i + wind_half; 
			
			if ( ((i % 1E6) == 0) )// || ((i > 61866584) && (i < 61868850)) )
			{
				System.err.println(target_chr + " " + i + "-" + to + "\tc=" + wind_c + " g=" + wind_g + " cg="+ wind_cg +
						" mid=" + middle + " exp_cg=" + exp_cg + " density=" + density);
			}
		
				if (middle >= 1)
				{
					// Adding the string version of the chromosome makes it quicker
					this.addScore(target_chr, middle+start_coord, new Double(density));
				}
		
		}
	}
	

//	// Returns a list of ChromFeaturePopulatorsRegion elements for the chromosome
//	public List wigFileRetrieveRegions(String wig_fn, String target_chr)
//	throws Exception
//	{
//		int chr_id = (new ChromFeatures()).chrom_from_public_str(target_chr); 
//			
//		// Returns a list of [chrom_id, start, end, strand, score, line_number] (float) objects
//		List regions = ChromFeaturePopulators.wigRegions(wig_fn, target_chr);
//		
//		// Iterate through them, replace the score with our score
//		Iterator it = regions.iterator();
//		while (it.hasNext())
//		{
//			ChromFeaturePopulatorsRegion reg = (ChromFeaturePopulatorsRegion)it.next();
//			
//			double cpg = ((Double)this.getScore(chr_id, reg.f_start)).doubleValue();
//			reg.f_score = cpg;
//		}
//	
//		return regions;
//	}
	
	

}
