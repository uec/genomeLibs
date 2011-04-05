package org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel;

public class BaseUtilsMore {

	static public byte[] toUpperCase(byte[] in)
	{
		byte[] out = new byte[in.length];
		for (int i = 0; i < in.length; i++)
		{
			out[i] = toUpperCase(in[i]);
		}
		return out;
	}

	static public byte toUpperCase(byte b)
	{
		return (byte)Character.toUpperCase((char)b);
	}
}
