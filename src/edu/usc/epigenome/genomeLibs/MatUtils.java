package edu.usc.epigenome.genomeLibs;

import java.util.*;
import java.io.*;

import org.apache.commons.math.stat.StatUtils;

public class MatUtils {

	
	public static int sumAll(int[][] mat)
	{
		return (int)(sumAll(MatUtils.intMatToDouble(mat)));
	}
	
	public static double sumAll(double[][] mat)
	{
	   	int nrow = mat.length;
    	int ncol = (mat[0]).length;

    	double out = 0.0;
    	for (int i = 0; i < nrow; i++)
    	{
    		for (int j = 0; j < ncol; j++)
    		{
    			out += mat[i][j];
    		}
    	}
    	
    	return out;
	}

	// Too slow to put it in a string
	public static void matlabCsv(PrintWriter writer, double[][] mat)
	{
		int n_r = mat.length;
		
		for (int i=0 ; i<n_r; i++)
		{
			int n_c = mat[i].length;
			for (int j=0; j<n_c; j++)
			{
				if (j>0) writer.print(",");
				writer.print(mat[i][j]);
				//System.err.println("mat["+i+"]["+j+"]="+ mat[i][j]);
			}
			writer.print("\n");
		}
	}
	
	
	public static String matString(int[][] mat)
	{
	   	int nrow = mat.length;
    	int ncol = (mat[0]).length;

    	String out = "";
    	for (int i = 0; i < nrow; i++)
    	{
    		if (i>0) out += "\n";
    		
    		for (int j = 0; j < ncol; j++)
    		{
    			if (j>0) out += "\t";
    			out += mat[i][j];
    		}
    	}
    	
    	return out;
	}	
	
	public static String matString(double[][] mat)
	{
	   	int nrow = mat.length;
    	int ncol = (mat[0]).length;

    	String out = "";
    	for (int i = 0; i < nrow; i++)
    	{
    		if (i>0) out += "\n";
    		
    		for (int j = 0; j < ncol; j++)
    		{
    			if (j>0) out += "\t";
    			out += mat[i][j];
    		}
    	}
    	
    	return out;
	}	

	public static double[] colMeans(double[][] m)
	{
		return colMeans(m,0);
	}

	public static double[] colMeans(double[][] m, int n_rows)
	{
		double[][] m_trans = MatUtils.transposedMat(m);
		return rowMeans(m_trans, n_rows);
	}
	
	public static double[] colSums(double[][] m)
	{
		double[][] m_trans = MatUtils.transposedMat(m);
		return rowSums(m_trans);
	}

	public static double[] colStdevs(double[][] m)
	{
		double[][] m_trans = MatUtils.transposedMat(m);
		return rowStdevs(m_trans);
	}

	public static double[] rowMeans(double[][] m)
	{
		return rowMeans(m,0);
	}
	
	public static double[] rowMeans(double[][] m, int n_cols)
	{
		int n_rows = m.length;
		double[] out = new double[n_rows];
		
		if (n_cols ==0) n_cols = m[0].length;
		for (int i = 0; i < n_rows; i++)
		{
			out[i] = StatUtils.mean(m[i], 0, n_cols);
		}
		return out;
	}
	
	public static double[] rowSums(double[][] m)
	{
		int n_rows = m.length;
		double[] out = new double[n_rows];
		
		for (int i = 0; i < n_rows; i++)
		{
			out[i] = StatUtils.sum(m[i]);
		}
		return out;
	}

	public static double[] rowStdevs(double[][] m)
	{
		int n_rows = m.length;
		double[] out = new double[n_rows];
		
		for (int i = 0; i < n_rows; i++)
		{
			out[i] = Math.sqrt(StatUtils.variance(m[i]));
		}
		return out;
	}
	
	protected static void downscaleArray(double[] out, double[] in)
	throws Exception
	{
		downscaleArray(out,in,0,out.length-1);
	}

	protected static void downscaleArray(double[] out, double[] in, int out_start_ind, int out_end_ind)
	throws Exception
	{
		int n_in = in.length;
		int n_out = out_end_ind - out_start_ind + 1;
//		if (n_in < n_out) System.err.println("downscaleArray(out,in) out must be smaller than in");

//		System.err.println("downscaleArray(" +n_out +"," + n_in + ", "+ out_start_ind + "-" + out_end_ind + ")");
		
		double factor = n_in / n_out;

		int next_in_pos = 0;
		for (int out_pos = 0; out_pos < n_out; out_pos++)
		{
			int out_ind = out_start_ind + out_pos;

			int in_s = Math.min(n_in-1,(int)(Math.round((out_pos) * factor))); // next_in_pos;
			int in_e = Math.min(n_in-1,(int)(Math.round((out_pos+1) * factor)));

			if (in_e >= in_s)
			{
				out[out_ind] = StatUtils.mean(in, in_s, in_e-in_s+1);
//				System.err.println("\tout[" + out_ind + "] = mean(" + in_s + "-" + in_e +")");
			}

//			next_in_pos = in_e+1;
		}

	}

	
	public static int[][] sumMats(int[][] m1, int[][] m2)
	{
		return MatUtils.doubleMatToInt(sumMats (MatUtils.intMatToDouble(m1), MatUtils.intMatToDouble(m2)));
	}
	
	public static double[][] sumMats(double[][] m1, double[][] m2)
	{
		int nrow = m1.length;
		int ncol = m1[0].length;
		
		double[][][] all_mats = new double[2][nrow][ncol];
		all_mats[0] = m1;
		all_mats[1] = m2;
		return sumMats(all_mats);
	}
	
	
	public static double[][] sumMats(double[][][] mats)
	{
		int nmats = mats.length;
		int nrow = mats[0].length;
		int ncol = mats[0][0].length;
		
		double[][] out = new double[nrow][ncol];
		for (int m = 0; m < nmats; m++)
		{
			for (int i = 0; i < nrow; i++)
			{
				for (int j = 0; j < ncol; j++)
				{
					out[i][j] += mats[m][i][j];
				}
			}
		}
		
		return out;
	}

	public static double[][] divMats(int[][] nmat , int[][] dmat)
	throws Exception
	{
		return divMats (MatUtils.intMatToDouble(nmat), MatUtils.intMatToDouble(dmat));
	}
		
	public static double[][] divMats(double[][] nmat , double[][] dmat)
	throws Exception
	{
		int nrow = nmat.length;
		int ncol = nmat[0].length;
		
		int drow = dmat.length;
		int dcol = dmat[0].length;

		if ((nrow != drow ) || (ncol != dcol))
		{
			throw new Exception("MatUtils::divMats(): nmat ("+ nrow + "," + ncol + ") doesn't have same dimensions as dmat ("+ drow + "," + dcol +")");
		}
		

		double[][] out = new double[nrow][ncol];
		for (int i = 0; i < nrow; i++)
		{
			for (int j = 0; j < ncol; j++)
			{
				out[i][j] = nmat[i][j] / dmat[i][j];
			}
		}
		
		return out;
	}

	
	public static void incrementMat(int[][] total , int[][] add)
	throws Exception
	{
		int nrow = total.length;
		int ncol = total[0].length;
		
		int drow = add.length;
		int dcol = add[0].length;

		if ((nrow != drow ) || (ncol != dcol))
		{
			throw new Exception("MatUtils::divMats(): nmat ("+ nrow + "," + ncol + ") doesn't have same dimensions as dmat ("+ drow + "," + dcol +")");
		}
		

		double[][] out = new double[nrow][ncol];
		for (int i = 0; i < nrow; i++)
		{
			for (int j = 0; j < ncol; j++)
			{
				total[i][j] += add[i][j];
			}
		}
	}

	public static int[][] invertedMat(int[][] in)
	{
		int nrow = in.length;
		int ncol = in[0].length;
		
		int[][] out = new int[nrow][ncol];
		for (int i = 0; i < nrow; i++)
		{
			for (int j = 0; j < ncol; j++)
			{
				out[(nrow-1)-i][(ncol-1)-j] = in[i][j];
			}
		}
		
		//System.err.println("Orig mat:\n" + MatUtils.matString(in) + "\nnew:\n" + MatUtils.matString(out) +"\n");
		
		return out;
	}
	
	public static double[][] transposedMat(double[][] in)
	{
		int nrow = in.length;
		int ncol = in[0].length;
		
		double[][] out = new double[ncol][nrow];
		for (int i = 0; i < nrow; i++)
		{
			for (int j = 0; j < ncol; j++)
			{
				out[j][i] = in[i][j];
			}
		}
		
	//	System.err.println("Orig mat:\n" + MatUtils.matString(in) + "\nnew:\n" + MatUtils.matString(out) +"\n");
		
		return out;
	}
	
	public static double[][] invertedMat(double[][] in)
	{
		int nrow = in.length;
		int ncol = in[0].length;
		
		double[][] out = new double[nrow][ncol];
		for (int i = 0; i < nrow; i++)
		{
			for (int j = 0; j < ncol; j++)
			{
				out[(nrow-1)-i][(ncol-1)-j] = in[i][j];
			}
		}
		
//		System.err.println("Orig mat:\n" + MatUtils.matString(in) + "\nnew:\n" + MatUtils.matString(out) +"\n");
		
		return out;
	}
	
	
	public static void incrementMat(double[][] total , double[][] add)
	throws Exception
	{
		int nrow = total.length;
		int ncol = total[0].length;
		
		int drow = add.length;
		int dcol = add[0].length;

		if ((nrow != drow ) || (ncol != dcol))
		{
			throw new Exception("MatUtils::divMats(): nmat ("+ nrow + "," + ncol + ") doesn't have same dimensions as dmat ("+ drow + "," + dcol +")");
		}
		

		double[][] out = new double[nrow][ncol];
		for (int i = 0; i < nrow; i++)
		{
			for (int j = 0; j < ncol; j++)
			{
				total[i][j] += add[i][j];
			}
		}
	}

	public static double[][] intMatToDouble(int[][] mat)
	{
	   	int nrow = mat.length;
    	int ncol = (mat[0]).length;

    	double[][] out = new double[nrow][ncol];
    	for (int i = 0; i < nrow; i++)
    	{
    		for (int j = 0; j < ncol; j++)
    		{
    			out[i][j] = (double)mat[i][j];
    		}
    	}
    	return out;
	}
	
	public static int[][] doubleMatToInt(double[][] mat)
	{
	   	int nrow = mat.length;
    	int ncol = (mat[0]).length;

    	int[][] out = new int[nrow][ncol];
    	for (int i = 0; i < nrow; i++)
    	{
    		for (int j = 0; j < ncol; j++)
    		{
    			out[i][j] = (int)mat[i][j];
    			//System.err.println(mat[i][j] + " -> " + out[i][j]);
    		}
    	}
    	return out;
	}
	
}
