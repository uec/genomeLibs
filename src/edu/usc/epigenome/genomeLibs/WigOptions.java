package edu.usc.epigenome.genomeLibs;

public class WigOptions {
	
	public String f_name = "generic wig";
	public int f_format = 0; // 0 (wig full), 1 (wig compact), 2 (GADA), 3 (BED)
	public boolean f_print_head = true; 
	public int f_step = 500;
	public int f_span = 500;
	public boolean f_upside_down = false;
	public int f_min_score = 0;
	String f_color = null;
	public int[] f_view_limits = null;
	
	public WigOptions()
	{
	}
	
	public WigOptions(String name, int format, int step, int span, 
			boolean upside_down, boolean print_head, String color,
			int min_score)
	{
		f_step = step;
		f_format = format;
		f_span = span;
		f_upside_down = upside_down;
		f_color = color;
		f_name = name;
		f_min_score = min_score;
	}

	public void makeFwStrand()
	{
		f_color = "0,0,255";
		f_upside_down = false;
	}
	
	public void makeRevStrand()
	{
		f_color = "255,0,0";
		f_upside_down = true;
	}
	
	public String trackHead()
	{
		String desc = f_name;
		String head = "track type=wiggle_0 name=\"" + f_name + "\" description=\"" + desc +
		"\" visibility=2 maxHeightPixels=128:25:11 graphType=bar windowingFunction=mean autoScale=off";
		
		if (f_view_limits != null)
		{
			int vl1 = (f_upside_down) ? (-1*f_view_limits[1]) : f_view_limits[0];
			int vl2 = (f_upside_down) ? (-1*f_view_limits[0]) : f_view_limits[1];
			head += " viewLimits=" + vl1 + ":" + vl2;
		}
		
		if (f_color != null)
		{
			head += " color=" + f_color;
		}
		
		return head;
	}
	
	public String fixedStepHead(String chrom, int chrom_start)
	{
		String out = "";
		out += "fixedStep\tchrom=" + chrom + "\tstart=" + chrom_start + 
		"\tstep=" + f_step;
		
		if (f_span > 0)
		{
			out += "\tspan=" + f_span;
		}

		return out;
	}
	
}
