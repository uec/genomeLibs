package edu.usc.epigenome.genomeLibs;

import java.io.*;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Logger;

import org.apache.commons.math.stat.StatUtils;

public class MatUtils {
	
	protected static Map<String,double[]> cachedExponentialWeights = new TreeMap<String,double[]>();

	public static double nanMean(double[] arr)
	{
		return nanMean(arr, 0, arr.length);
	}
	
	public static double nanMin(double[] arr)
	{
		double val = Double.NaN;
		for (int i = 0; i < arr.length; i++)
		{
			if (Double.isNaN(val))
			{
				val = arr[i];
			}
			else if (Double.isNaN(arr[i]))
			{
				// Do nothing	
			}
			else
			{
				if (arr[i]<val) val = arr[i];
			}
		}
		
		return val;
	}

	public static double nanMax(double[] arr)
	{
		double val = Double.NaN;
		for (int i = 0; i < arr.length; i++)
		{
			if (Double.isNaN(val))
			{
				val = arr[i];
			}
			else if (Double.isNaN(arr[i]))
			{
				// Do nothing	
			}
			else
			{
				if (arr[i]>val) val = arr[i];
			}
		}
		
		return val;
		
	}

	
	public static double nanMean(double[] arr, int startInd, int len)
	{
		double total = 0.0;
		double count = 0.0;
		
		for (int i = startInd; i < (startInd+len) ; i++)
		{
			double val = arr[i];
			//System.err.printf("\t\tFound val=%f\n",val);
			if (!Double.isNaN(val))
			{
				total += val;
				count += 1.0;
			}
		}
		
		return (count == 0.0) ? Double.NaN : (total/count);
	}
	
	// exp=-1/3 and distMin=10 are good
	public static double nanMeanExponentialWeighting(double[] arr, double exponentialFactor, int distMin)
	{
		int startInd = 0;
		int len = arr.length;
		double total = 0.0;
		double count = 0.0;
		
		// Cache the weights so we don't have to calculate a million times
		String key = String.format("%d__%.4f",len,exponentialFactor);
		double[] weights = MatUtils.cachedExponentialWeights.get(key);
		if (weights==null)
		{
			Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Caching exponential weights for " + key);
			double mid = Math.round((double)len/2.0);
			double maxWeight = Math.pow(distMin, exponentialFactor);
			weights = new double[len];
			for (int i = 0; i < len; i++)
			{
				double dist = Math.abs((double)i - mid);
				if (dist<(double)distMin) dist = (double)distMin;
				weights[i] = Math.pow(dist, exponentialFactor) / maxWeight;
				Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine(String.format("weights[%d]=%.4f",i,weights[i]));
			}
			
			MatUtils.cachedExponentialWeights.put(key, weights);
		}
		
		for (int i = startInd; i < (startInd+len) ; i++)
		{
			double val = arr[i];
			//System.err.printf("\t\tFound val=%f\n",val);
			if (!Double.isNaN(val))
			{
				total += (val*weights[i]);
				count += 1.0;
			}
		}
		
		return (count == 0.0) ? Double.NaN : (total/count);
	}
	
	public static double nanStdev(double[] arr)
	{
		return nanStdev(arr, 0, arr.length);
	}

	public static double nanStdev(double[] arr , int startInd, int len)
	{
		return Math.sqrt(nanVariance(arr, startInd, len));
	}
	
	
	public static double nanVariance(double[] arr)
	{
		return nanVariance(arr, 0, arr.length);
	}
	
	
	public static double nanVariance(double[] arr , int startInd, int len)
	{
		double mean = nanMean(arr, startInd, len);
		
		int count = 0;
		double totalVar = 0;
		for (int i = startInd; i < (startInd+len) ; i++)
		{
			if (!Double.isNaN(arr[i]))
			{
				count++;
				totalVar += Math.pow(arr[i] - mean, 2.0);
			}
		}		
		
		return totalVar / (double)count;
	}
	
	
	public static double nanSum(double[] arr)
	{
		return nanSum(arr, 0, arr.length);
	}

	public static double nanSum(double[] arr, int startInd, int len)
	{
		double total = 0.0;
		boolean found = false;
		
		for (int i = startInd; i < (startInd+len) ; i++)
		{
			double val = arr[i];
			if (!Double.isNaN(val))
			{
				total += val;
				found = true;
			}
		}
		
		return (!found) ? Double.NaN : total;
	}


		
	public static void initMat(int[] mat, int initVal)
	{
		for (int i = 0; i < mat.length; i++)
		{
			mat[i] = initVal;
		}
	}


	public static void nansToVal(double[] mat, double inVal)
	{
		for (int i = 0; i < mat.length; i++)
		{
			if (Double.isNaN(mat[i])) mat[i] = inVal;
		}
	}

	public static int countNans(double[][] mat)
	{
		int total = 0;
		for (int i = 0; i < mat.length; i++)
		{
			for (int j = 0; j < mat[0].length; j++)
			{
				if (Double.isNaN(mat[i][j])) total++;
			}
		}
		return total;
	}
	
	public static void nansToVal(double[][] mat, double inVal)
	{
		for (int i = 0; i < mat.length; i++)
		{
			for (int j = 0; j < mat[0].length; j++)
			{
				if (Double.isNaN(mat[i][j])) mat[i][j] = inVal;
			}
		}
	}

	public static void initMat(double[] mat, double initVal)
	{
		for (int i = 0; i < mat.length; i++)
		{
			mat[i] = initVal;
		}
	}
	
	public static void initMat(int[][] mat, int initVal)
	{
		for (int i = 0; i < mat.length; i++)
		{
			for (int j = 0; j < mat[0].length; j++)
			{
				mat[i][j] = initVal;
			}
		}
	}


	public static void initMat(double[][] mat, double initVal)
	{
		for (int i = 0; i < mat.length; i++)
		{
			for (int j = 0; j < mat[0].length; j++)
			{
				mat[i][j] = initVal;
			}
		}
	}

	public static void initMat(int[][][] mat, int initVal)
	{
		for (int i = 0; i < mat.length; i++)
		{
			for (int j = 0; j < mat[0].length; j++)
			{
				for (int k = 0; k < mat[0][0].length; k++)
				{
					mat[i][j][k] = initVal;
				}
			}
		}
	}
	
	public static void initMat(double[][][] mat, double initVal)
	{
		for (int i = 0; i < mat.length; i++)
		{
			for (int j = 0; j < mat[0].length; j++)
			{
				for (int k = 0; k < mat[0][0].length; k++)
				{
					mat[i][j][k] = initVal;
				}
			}
		}
	}
	
	
	// exponetial weighting goes from the mid point out.
	public static double[][] sortRows(double[][] in)
	{
		return sortRows(in,0.0,0);
	}
	
	public static double[][] sortRows(double[][] in, double exponentialFactor, int minDist)
	{
		int nr = in.length;
		int nc = in[0].length;
		
		double[][] out = new double[nr][];
		SortedMap<Double,double[]> sorter = new TreeMap<Double,double[]>();
		for (int i = 0; i < nr; i++)
		{
			double mean = (exponentialFactor!=0.0) ?  
					MatUtils.nanMeanExponentialWeighting(in[i], exponentialFactor,minDist) : MatUtils.nanMean(in[i]);
			
			// A map will lose rows with identical means.  So we have to do some
			// tricks to keep keys unique
			Double key = (Double.isNaN(mean)) ? new Double(-1 * 1E-5) : new Double(mean);
			while (sorter.containsKey(key))
			{
				Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).fine("Got duplicate key " + key + "... incrementing\n");
				key += 1E-8;
			}
			
			sorter.put(new Double(key),in[i]);
		}
		
		int i = 0;
		for (double[] row : sorter.values()) // Guaranteed to be sorted
		{
			out[i] = row;
			i++;
		}
		
		// Double check that we're not making a ragged matrix
		if (i != nr)
		{
			System.err.printf("MatUtils.sortRows got a matrix with %d rows, but sorter only returned %d\n", nr,i);
			System.exit(1);
		}
		
		return out;
	}
	

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


	
	public static void matlabCsv(PrintWriter writer, double[][] mat)
	{
		matlabCsv(writer, mat, 0, 0);
	}
		
	// If nI or nJ are 0, it uses the total number of rows/columns in the array
	public static void matlabCsv(PrintWriter writer, double[][] mat, int nI, int nJ)
	{
		int n_r = mat.length;
		if (nI>0 && nI<=n_r) n_r = nI;
		int n_c = mat[0].length;
		if (nJ>0 && nJ<=n_c) n_c = nJ;
		
		for (int i=0 ; i<n_r; i++)
		{
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
			out[i] = nanMean(m[i], 0, n_cols);
		}
		return out;
	}
	
	public static double[] vectSum(double[] a, double[] b)
	{
		int nC = a.length;
		double[][] rowMat = new double[2][nC];
		rowMat[0] = a;
		rowMat[1] = b;
		return MatUtils.colSums(rowMat);
	}
	
	public static double[] colSums(double[][] m)
	{
		double[][] m_trans = MatUtils.transposedMat(m);
		return rowSums(m_trans);
	}

	public static double[] rowSums(double[][] m)
	{
		int n_rows = m.length;
		double[] out = new double[n_rows];
		
		for (int i = 0; i < n_rows; i++)
		{
			out[i] = nanSum(m[i]);
		}
		return out;
	}

	public static double[] rowStdevs(double[][] m)
	{
		int n_rows = m.length;
		double[] out = new double[n_rows];
		
		for (int i = 0; i < n_rows; i++)
		{
			out[i] = Math.sqrt(nanVariance(m[i]));
		}
		return out;
	}
	
	public static double[][] downscaleMat(double[][] mat, int numRowsNew, int numColsNew, double smoothingFact)
	throws Exception
	{
		int n = mat.length;
		int m = mat[0].length;
		
		// This isn't very efficient, especially memory wise.  But it's simple code.  
		// First rows
		if ((numColsNew!=0) && (m!=numColsNew)) mat = MatUtils.downscaleMatRows(mat, numColsNew,  smoothingFact);
		
		// Then cols
		if  ((numRowsNew!=0) && (n!=numRowsNew)) mat = MatUtils.downscaleMatCols(mat, numRowsNew,  smoothingFact);
		
		return mat;
	}
	
	public static double[][] downscaleMatCols(double[][] in, int numRowsNew, double smoothingFact)
	throws Exception
	{
		// Probably not so efficient , but hey.
		double[][] out = MatUtils.transposedMat(MatUtils.downscaleMatRows(MatUtils.transposedMat(in), numRowsNew, smoothingFact));
		return out;
	}
	
	public static double[][] downscaleMatRows(double[][] in, int numColsNew, double smoothingFact)
	throws Exception
	{
		int n = in.length;
		int m = in[0].length;
		
		int m1 = numColsNew;
		
		if (m == m1) return in;
		
		double[][] out = new double[n][];
		for (int i = 0; i < n; i++)
		{
			double[] newr = new double[m1];
			downscaleArray(newr, in[i], smoothingFact);
			out[i] = newr;
		}
		
		return out;
	}
	
	
	public static void downscaleArray(double[] out, double[] in, double smoothingFact)
	throws Exception
	{
		downscaleArray(out,in,0,out.length-1, smoothingFact);
	}

	public static void downscaleArray(double[] out, double[] in)
	throws Exception
	{
		downscaleArray(out,in,0,out.length-1);
	}

	public static void downscaleArray(double[] out, double[] in, int out_start_ind, int out_end_ind)
	throws Exception
	{
		downscaleArray(out,in,out_start_ind,out_end_ind,0.0);
	}
	
	public static void downscaleArray(double[] out, double[] in, int out_start_ind, int out_end_ind, double SMOOTHING_FACT)
	throws Exception
	{
		
		int n_in = in.length;
		int n_out = out_end_ind - out_start_ind + 1;
//		if (n_in < n_out) System.err.println("downscaleArray(out,in) out must be smaller than in");

//		System.err.println("downscaleArray(" +n_out +"," + n_in + ", "+ out_start_ind + "-" + out_end_ind + ")");
		
		double factor = n_in / n_out;
		
		double smoothingHalf = Math.ceil(SMOOTHING_FACT/2);

//		int next_in_pos = 0;
		for (int out_pos = 0; out_pos < n_out; out_pos++)
		{
			int out_ind = out_start_ind + out_pos;

			int in_s = Math.min(n_in-1,(int)(Math.round((out_pos - smoothingHalf) * factor))); // next_in_pos;
			int in_e = Math.min(n_in-1,(int)(Math.round((out_pos+1+smoothingHalf) * factor)));

			in_s = Math.max(0, in_s);
			in_e = Math.min(n_in-1,in_e);
			
			if (in_e >= in_s)
			{
				out[out_ind] = nanMean(in, in_s, in_e-in_s+1);
				//System.err.printf("\tout[%d]=mean(in[%d-%d])=%f\n",out_ind,in_s,in_e,out[out_ind]);
			}

//			next_in_pos = in_e+1;
		}

	}

	
	public static int[][] nanMeanMats(int[][] m1, int[][] m2)
	{
		return MatUtils.doubleMatToInt(nanMeanMats (MatUtils.intMatToDouble(m1), MatUtils.intMatToDouble(m2)));
	}
	
	public static double[][] nanMeanMats(double[][] m1, double[][] m2)
	{
		int nrow = m1.length;
		int ncol = m1[0].length;
		
		double[][][] all_mats = new double[2][nrow][ncol];
		all_mats[0] = m1;
		all_mats[1] = m2;
		return nanMeanMats(all_mats);
	}
	
	
	public static double[][] nanMeanMats(double[][][] mats)
	{
		int nmats = mats.length;
		int nrow = mats[0].length;
		int ncol = mats[0][0].length;

		double[][] out = new double[nrow][ncol];

		// We can re-use this.
		double[] series = new double[nmats];

		
		for (int i = 0; i < nrow; i++)
		{
			for (int j = 0; j < ncol; j++)
			{
				// For nan mean, the denominator will be the number of non-Nans
				for (int m = 0; m < nmats; m++)
				{
					series[m] = mats[m][i][j];
				}
				out[i][j] = MatUtils.nanMean(series);
			}
		}

		return out;
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
		
		// If input mats have NaNs, we are conservative
		// and initialize to nan
		MatUtils.initMat(out, Double.NaN);


		for (int m = 0; m < nmats; m++)
		{
			for (int i = 0; i < nrow; i++)
			{
				for (int j = 0; j < ncol; j++)
				{

					if (Double.isNaN(out[i][j]))
					{
						out[i][j] = mats[m][i][j];
					}
					else
					{
						if (!Double.isNaN(mats[m][i][j])) out[i][j] += mats[m][i][j];
					}
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
	
	public static double[] divVects(double[] a, double[] b)
	throws Exception
	{
		double[][] aMat = new double[1][a.length];
		aMat[0] = a;
		double[][] bMat = new double[1][b.length];
		bMat[0] = b;
		
		double[][] outMat = MatUtils.divMats(aMat, bMat);
		return outMat[0];
	}
	
	public static double[] divVect(double[] vect, double denom)
	{
		return multVect(vect, 1/denom);
	}

	public static double[] multVect(double[] vect, double denom)
	{
		
		double[][] mat = new double[1][vect.length];
		mat[0] = vect;
		double[][] outMat = multMat(mat, denom);
		return outMat[0];
	}
	
	public static double[][] divMat(double[][] mat, double denom)
	{
		return multMat(mat, 1/denom);
	}

	public static double[][] multMat(double[][] mat, double denom)
	{
		int nrow = mat.length;
		int ncol = mat[0].length;

		double[][] out = new double[nrow][ncol];
		for (int i = 0; i < nrow; i++)
		{
			for (int j = 0; j < ncol; j++)
			{
				out[i][j] = mat[i][j] * denom;
			}
		}
		
		return out;		
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
		
		System.err.printf("Transposing mat(%d,%d)\n",nrow,ncol);
		
		double[][] out = new double[ncol][nrow];
		for (int i = 0; i < nrow; i++)
		{
			for (int j = 0; j < ncol; j++)
			{
				try
				{
					out[j][i] = in[i][j];
				}
				catch (Exception e)
				{
					System.err.printf("out[%d][%d]=in[%d][%d]\n(in[%d])=%s\nArray probably ragged\n",
							j,i,i,j,i,""+in[i]);
					System.exit(1);
				}
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
		

		for (int i = 0; i < nrow; i++)
		{
			for (int j = 0; j < ncol; j++)
			{
				if (Double.isNaN(total[i][j]))
				{
					total[i][j] = add[i][j];
				}
				else
				{
					if (!Double.isNaN(add[i][j])) total[i][j] += add[i][j];
				}
			}
		}
	}

	public static double[][][] intMatToDouble(int[][][] mat)
	{
		double[][][] out = new double[mat.length][][];
		for (int i = 0; i < mat.length; i++)
		{
			out[i] = MatUtils.intMatToDouble(mat[i]);
		}
		return out;
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
