/**
 * 
 */
package edu.usc.epigenome.genomeLibs.Counters;

import java.util.HashMap;

import org.biojava.bio.symbol.Symbol;

import edu.usc.epigenome.genomeLibs.ReadPos.ReadPos;
import edu.usc.epigenome.genomeLibs.ReadPos.ReadPosRich;


/**
 * A collection of counts for a ReadPos
 * @author benb
 *  12-09-08 zack: added google charting and keyset selectors
 *
 */


public class ReadPosCounter extends TreeMapCounter<ReadPos>  {


	private static final long serialVersionUID = 898341971045049360L;

	/**
	 * 
	 */
	public ReadPosCounter() {
	}
	
	 
	/**
	 * google bar chart of counts by symbol (all cycles totaled)
	 * @return an url of the chart image
	 */
	public String getSymbolCountChartURL()
	{
		String dataString = "&chd=t:";
		String labelString = "&chl=";
		HashMap<String,Integer> counts = new HashMap<String,Integer>();
		long min = 0;
		long max = 0;
		for(ReadPos key :  this.keySet())
		{
			if(counts.containsKey(String.valueOf(key.getSymToken())))
				counts.put(String.valueOf(key.getSymToken()), counts.get(String.valueOf(key.getSymToken())) + this.get(key));
			else
				counts.put(String.valueOf(key.getSymToken()), this.get(key));
		}
		
		for(String Skey :  counts.keySet())
		{
			labelString += Skey + "|";
			dataString += counts.get(Skey) + ",";
			if(counts.get(Skey) > max)
				max = counts.get(Skey);
			if(counts.get(Skey) < min)
				min = counts.get(Skey);
		}
		
		String retVal = "http://chart.apis.google.com/chart?chs=530x520&chco=ff0000|00ff00|0000ff|ff00ff|00ffff&chxt=y&cht=bvs&chxr=0," + min +"," + max + "&chds=" + min +"," + max + dataString.substring(0, dataString.length() - 1) + labelString; 
		return retVal;
	}
	
	
	/**
	 * google chart of symbols by count per cycle
	 * @return an url of the chart image
	 */
	public String getCountByCyclesChartURL()
	{
		long min = 0;
		long max = 0;
		long cyclesMax=0;
		String dataString = "&chd=t:";
		String labelString = "&chdl=";
		String scalingString = "&chds=";
		String lineStyleString = "&chls=";
		HashMap<String,String> countsDataGroup = new HashMap<String,String>();
		HashMap<String,String> cyclesDataGroup = new HashMap<String,String>();
		
		for(ReadPos tkey : this.keySet())
		{
			ReadPosRich key = (ReadPosRich) tkey;
			if(this.get(key) > max)
				max = this.get(key);
			if(this.get(key) < min)
				min = this.get(key);
			if(key.getCycle() > cyclesMax)
				cyclesMax = key.getCycle();
			

			if(countsDataGroup.containsKey(String.valueOf(key.getSymToken())))
			{
				countsDataGroup.put(String.valueOf(key.getSymToken()), countsDataGroup.get(String.valueOf(key.getSymToken())) + "," + this.get(key));
				cyclesDataGroup.put(String.valueOf(key.getSymToken()), cyclesDataGroup.get(String.valueOf(key.getSymToken())) + "," + key.getCycle());
			}
			else
			{
				countsDataGroup.put(String.valueOf(key.getSymToken()), String.valueOf(this.get(key)));
				cyclesDataGroup.put(String.valueOf(key.getSymToken()),  String.valueOf(key.getCycle()));
			}						
		}
		for(String datakey : countsDataGroup.keySet())
		{
			labelString += datakey + "|";
			//dataString += counts.get(key) + ",";
			dataString += cyclesDataGroup.get(datakey) + "|" + countsDataGroup.get(datakey) + "|";
			scalingString += "0," + cyclesMax + "," + min + "," + max + ",";
			lineStyleString += "2,0,0|";					
		}
		double gridLine = 100.0 / cyclesMax;
		dataString = dataString.substring(0, dataString.length() - 1);
		labelString = labelString.substring(0,labelString.length() - 1);
		scalingString = scalingString.substring(0,scalingString.length() -1);
		lineStyleString = lineStyleString.substring(0,lineStyleString.length() -1);
		String retVal = "http://chart.apis.google.com/chart?&cht=lxy&chs=530x530&chco=ff0000,00dd00,0000ff,dd00dd,00dddd,dddd00&chg=" + gridLine +",0.0,2.0,1&chxt=x,y&chxr=0," + 0 +"," + cyclesMax + "|1," + min +"," + max + lineStyleString + scalingString + dataString + labelString;
		return retVal;		
	}
	

	/**
	 * generate a google chart of the symbols by cycle, normalized by total counts per cycle
	 * @return a google chart image URL
	 */
	public String getCountByCyclesPercentageChartURL()
	{
		long min = 0;
		long max = 100;
		long cyclesMax=0;
		String dataString = "&chd=t:";
		String labelString = "&chdl=";
		String scalingString = "&chds=";
		String lineStyleString = "&chls=";
		HashMap<String,String> countsDataGroup = new HashMap<String,String>();
		HashMap<String,String> cyclesDataGroup = new HashMap<String,String>();
		
		for(ReadPos tkey : this.keySet())
		{
			ReadPosRich key = (ReadPosRich) tkey;
			if(key.getCycle() > cyclesMax)
				cyclesMax = key.getCycle();
			long total = 0;
			ReadPosCounter byCycleCount = this.getKeysByCycle(key.getCycle());
			for(ReadPos cycleKey : byCycleCount.keySet())
			{
				total += byCycleCount.get(cycleKey);
			}

			if(countsDataGroup.containsKey(String.valueOf(key.getSymToken())))
			{
				countsDataGroup.put(String.valueOf(key.getSymToken()), countsDataGroup.get(String.valueOf(key.getSymToken())) + "," +  (100L * this.get(key) / total));
				cyclesDataGroup.put(String.valueOf(key.getSymToken()), cyclesDataGroup.get(String.valueOf(key.getSymToken())) + "," + key.getCycle());
			}
			else
			{
				countsDataGroup.put(String.valueOf(key.getSymToken()), String.valueOf(( 100L * this.get(key) / total)));
				cyclesDataGroup.put(String.valueOf(key.getSymToken()),  String.valueOf(key.getCycle()));
			}						
		}
		for(String datakey : countsDataGroup.keySet())
		{
			labelString += datakey + "|";
			//dataString += counts.get(key) + ",";
			dataString += cyclesDataGroup.get(datakey) + "|" + countsDataGroup.get(datakey) + "|";
			scalingString += "0," + cyclesMax + "," + min + "," + max + ",";
			lineStyleString += "2,0,0|";					
		}
		double gridLine = 100.0 / cyclesMax;
		dataString = dataString.substring(0, dataString.length() - 1);
		labelString = labelString.substring(0,labelString.length() - 1);
		scalingString = scalingString.substring(0,scalingString.length() -1);
		lineStyleString = lineStyleString.substring(0,lineStyleString.length() -1);
		String retVal = "http://chart.apis.google.com/chart?&cht=lxy&chs=530x530&chco=ff0000,00dd00,0000ff,dd00dd,00dddd,dddd00&chg=" + gridLine +",0.0,2.0,1&chxt=x,y&chxr=0," + 0 +"," + cyclesMax + "|1," + min +"," + max + lineStyleString + scalingString + dataString + labelString;
		return retVal;		
	}
	
	
	/**
	 * subset of this collection that only pertains to the given cycle
	 * @param cycle only keys for this cycle value will be returned
	 * @return the new ReadPosCounter filtered
	 */
	protected ReadPosCounter getKeysByCycle(int cycle)
	{
		ReadPosCounter ret = new ReadPosCounter();
		for(ReadPos tkey : this.keySet())
		{
			ReadPosRich key = (ReadPosRich) tkey;
			if(key.getCycle() == cycle)
				ret.put(key, this.get(key));
		}
		return ret;
	}	
	
	/**
	 * subset of this collection that only pertains to the given symbol
	 * @param a symbol that will be used to filter all entries
	 * @return the new ReadPosCounter filtered
	 */
	protected ReadPosCounter getKeysBySymbol(Symbol symbol)
	{
		ReadPosCounter ret = new ReadPosCounter();
		for(ReadPos key : this.keySet())
		{
			if(key.getSym() == symbol)
				ret.put(key, this.get(key));
		}
		return ret;
	}
	
	 
	/**
	 * //get subset of this collection that has quality greater then specified
	 * @param quality value for the cutoff
	 * @return the new ReadPosCounter filtered
	 */
	protected ReadPosCounter getKeysGreaterThenQuality(int quality)
	{
		ReadPosCounter ret = new ReadPosCounter();
		for(ReadPos tkey : this.keySet())
		{
			ReadPosRich key = (ReadPosRich) tkey;
			if(key.getQual() > quality)
				ret.put(key, this.get(key));
		}
		return ret;
	}
	
	 
	/**
	 * get subset of this collection that has quality less then specified
	 * @param quality value for the cutoff
	 * @return the new ReadPosCounter filtered
	 */
	protected ReadPosCounter getKeysLessThenQuality(int quality)
	{
		ReadPosCounter ret = new ReadPosCounter();
		for(ReadPos tkey : this.keySet())
		{
			ReadPosRich key = (ReadPosRich) tkey;
			if(key.getQual() < quality)
				ret.put(key, this.get(key));
		}
		return ret;
	}
	
	
	/**
	 * Reduce and remove all cycle quality strand information and just
	 * return raw symbol counts
	 * @return
	 */
	public SymbolCounter getSymbolCounts()
	{
		SymbolCounter out = new SymbolCounter();
		for(ReadPos tkey : this.keySet())
		{
			Symbol sym = tkey.getSym();
			out.increment(sym);
		}

		return out;
	}

	
	
}
