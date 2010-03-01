package edu.usc.epigenome.genomeLibs;

import java.util.Map;
import java.util.TreeMap;

import com.mallardsoft.tuple.Pair;
import com.mallardsoft.tuple.Tuple;

public class DoubleRange extends Pair<Double,Double> {

	public DoubleRange(Double m1, Double m2) {
		super(m1, m2);
	}
	


}


//public double a;
//public double b;
//
///**
// * 
// */
//public DoubleRange(double inA, double inB) {
//	a = inA;
//	b = inB;
//}
//
//
//
//
//public static String rangeToString(double a, double b)
//{
//	DoubleRange dr = new DoubleRange(a,b);
//	return rangeToString(dr);
//}
//
//public static String rangeToString(DoubleRange dr)
//{
//	return dr.toString();
//}
//
//public static String rangePairToString(DoubleRange dra, DoubleRange drb)
//{
//	String out = rangeToStringMap.get(dr);
//	if (out == null)
//	{
//		out = dr.toString();
//		rangeToStringMap.put(dr, out);
//	}
//	return out;
//}
//
//
//
///* (non-Javadoc)
// * @see java.lang.Object#toString()
// */
//@Override
//public String toString() {
//	String out = rangeToStringMap.get(this);
//	if (out == null)
//	{
//		out = String.format("%f__%f", this.a, this.b);
//		rangeToStringMap.put(this, out);
//	}
//	return out;
//}
//
//
//
//
///* (non-Javadoc)
// * @see java.lang.Object#hashCode()
// */
//@Override
//public int hashCode() {
//	String str = this.toString();
//	return str.hashCode();		
//}
//
//
//
//
///* (non-Javadoc)
// * @see java.lang.Object#equals(java.lang.Object)
// */
//@Override
//public boolean equals(Object obj)
//{
//
//	boolean out = false;
//	try
//	{
//		DoubleRange other = (DoubleRange)obj;
//		if ((other.a != this.a) || (other.b != this.b)) out = false;
//	}
//	catch (Exception e)
//	{
//		System.err.println("Why are we comparing a DoubleRange to a " + obj.getClass());
//	}
//
//	return out;
//}

