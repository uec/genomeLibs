package edu.usc.epigenome.genomeLibs;

public class MiscUtils {


	public static String revCompNucStr(String in)
	{
		String out = "";
		
		for (int i = (in.length()-1); i >= 0; i--)
		{
			out += revCompNuc(in.charAt(i));
		}

		System.err.println("revCompStr(" + in + ") = " + out);
		return out;
	}

	public static char revCompNuc(char c)
	{
		char cout = c;

		switch (c) {
		case 'A': cout = 'T'; break;
		case 'a': cout = 't'; break;
		case 'T': cout = 'A'; break;
		case 't': cout = 'a'; break;
		case 'C': cout = 'G'; break;
		case 'c': cout = 'g'; break;
		case 'G': cout = 'C'; break;
		case 'g': cout = 'c'; break;
		}

		return cout;
	}
	
	public static char toggleCase(char c)
	{
		
		char out;
		if (Character.isUpperCase(c))
		{
			out = Character.toLowerCase(c);
		}
		else if (Character.isLowerCase(c))
		{
			out = Character.toUpperCase(c);
		}
		else
		{
			out = c;
		}
		return out;
	}
	
	
	public static String twodArrayString(int[][] a)
	{
		int ni = a.length;
		int nj = (a[0].length);
		
		String out = "";
		for (int i = 0; i < ni; i++)
		{
			for (int j = 0; j < nj; j++)
			{
				if (j>0) out += "\t";
				out += a[i][j];
			}
			out += "\n";
		}
		
		return out;
	}
	
//	// arg is a command line arg, which can be either a class name
//	// or a GFF file.  An unpopulated ChromFeatures object is 
//	// returned.  You must call populate() to populate the class.
//	public static ChromFeatures fetchChromFeaturesArg(String arg)
//	{
//
//		ChromFeatures feats = null;
//
//		String target_class_name = arg;
//		if (arg.indexOf(".") < 0)
//		{
//			target_class_name = "org.usckeck.genome." + target_class_name;
//		}
//
//		try
//		{
//			Class target_class = Class.forName(target_class_name);
//			feats = (ChromFeatures)target_class.newInstance();
//		}
//		catch (Exception e)
//		{
//			// Not a class
//		}
//
//		// If we're here , it's a file
//		if (feats == null)
//		{
//			feats= new ChromFeatures(arg);
//		}
//
//		return feats;
//	}
//	
//	
//	public static BitSet gtfBitmask(String fn, int chr_int)
//	throws Exception
//	{
//			ChromFeatures cf = new ChromFeatures(fn, false);
//			cf.c_first_chrom_num = chr_int;
//			cf.c_last_chrom_num = chr_int;
//			System.err.println("\tPopulating ChromFeatures " +fn);
//			cf.populate();
//			System.err.println("\tPopulated ChromFeatures, " + fn + " , " + cf.num_features() + " features");
//			return cf.bitSet(chr_int , false);
//	}
//
//	// filter_repmsk should be reponly or norep.  If it's anything else, we
//	// return null.
//	public static BitSet repeatBitmask(String filter_repmsk, int chr_int)
//	throws Exception
//	{
//		BitSet repmsk = null;
//		if (filter_repmsk.equalsIgnoreCase("reponly") || filter_repmsk.equalsIgnoreCase("norep"))
//		{
//			ChromFeatures repmat = new DbRepeatMasker();
////			ChromFeatures repmat = new DbRepeatMaskerLine();
//			repmat.c_first_chrom_num = chr_int;
//			repmat.c_last_chrom_num = chr_int;
//			System.err.println("\tPopulating DbRepeatMasker");
//			repmat.populate();
//			System.err.println("\tPopulated DbRepeatMasker, " + repmat.num_features() + " features");
//
//			ChromFeatures simple = null;
//			simple = new DbSimpleRepeat();
//			simple.c_first_chrom_num = chr_int;
//			simple.c_last_chrom_num = chr_int;
//			System.err.println("\tPopulating DbSimpleRepeat");
//			simple.populate();
//			System.err.println("\tPopulated DbSimpleRepeat, " + simple.num_features() + " features");
//
//			repmsk = repmat.bitSet(chr_int , false);
//			if (simple != null) repmsk.or(simple.bitSet(chr_int, false));
//			if (filter_repmsk.equalsIgnoreCase("norep")) repmsk.flip(0,repmsk.length()-1);
//		}
//
//		return repmsk;
//	}


	
	
}
