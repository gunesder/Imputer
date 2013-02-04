package edu.berkeley.path.imputer;

import org.apache.commons.math.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math.analysis.polynomials.PolynomialSplineFunction;

import scala.actors.threadpool.Arrays;

import java.lang.*;

public class MyUtilities {
	
	// *********************************************************************************************
	// ************ Array operations ***************************************************************
	// *********************************************************************************************
	
	// Scale Vector
	public static double[] scaleVector(double[] d, double s){
		for (int k=0;k<d.length;k++){
			d[k]*=s;			
		}
		return d;
	}
	
	// Scale Matrix
	public static double[][] scaleMatrix(double[][] M, double s){
		double[][] out = new double[M.length][M[0].length];
		for (int k=0;k<M.length;k++){
			out[k] = scaleVector(M[k],s);			
		}
		return out;
	}
	
	// Adding
	public static double[] addVectors(double[] a, double[] b){
		if (a.length != b.length){
			throw new ArithmeticException();
		}
		double[] sum = new double[a.length];
		for (int k=0;k<a.length;k++){
			sum[k] = a[k] + b[k];			
		}
		return sum;
		
	}
	
	// Assign Vector to a specific column of a 2D array
	public static double[][] assignColumn(double[][] matrix, double[] vector, int columnIndex){
		if (vector.length != matrix.length){
			throw new ArithmeticException();
		}
		for (int k=0;k<vector.length;k++){
			matrix[k][columnIndex] = vector[k];			
		}
		return matrix;
	}
	
	// Append Column Vector to 2D array
	public static double[][] appendColumn(double[][] matrix, double[] vector){
		if (vector.length != matrix.length){
			throw new ArithmeticException();
		}
		double[][] matrixOut = new double[matrix.length][matrix[0].length+1];
		for (int i = 0; i < matrixOut.length; i++) {
		    System.arraycopy(matrix[i], 0, matrixOut[i], 0, matrix[0].length);
		}
		for (int k=0;k<vector.length;k++){
			matrixOut[k][matrix[0].length] = vector[k];			
		}
		return matrixOut;
	}
	
	// Add columns of 2D array
	public static double[] sumColumns(double[][] matrix){
		double[] out = new double[matrix.length];
		for (int k=0;k<out.length;k++){
			for (int l=0;l<matrix[k].length;l++){
				out[k]+= matrix[k][l];				
			}			 
		}
		return out;
	}
	
	// Interpolate Vector
	public static double[] interpolateVector(double[] vector, int newsize){
		// allocate new vector with new size
		double[] resultVector = new double[newsize];
		// create array of interpolation points
		double[] interpolationPoints = MyUtilities.createIncrementVector(0,(double) vector.length-1,1);
		// interpolate and assign to output
		LinearInterpolator interliPolat = new LinearInterpolator();
		PolynomialSplineFunction resultFunctionSet = interliPolat.interpolate(interpolationPoints, vector);
		// evaluate the resulting interpolation function at the interpolationPoints
		PolynomialFunction[] functions = resultFunctionSet.getPolynomials();
		double[] knots = resultFunctionSet.getKnots();
		double[] evaluationPoints = MyUtilities.createIncrementVector(0,(double) vector.length-1,(double) vector.length/(newsize+1));
		for (int k = 0;k<evaluationPoints.length;k++){
			// determine which function to use (i.e. in which knot interval is the current evaluation point)
			int m;
			for (m = 0;m<knots.length;m++){
				if (evaluationPoints[k]<knots[m]){
					m--;
					break;
				} 
			}
			if (m == knots.length) m-=2;
			resultVector[k] = functions[m].value(evaluationPoints[k]-knots[m]); 
		}
		return resultVector;
	}
	
	// create m to n array (matlab 1:n)
	public static double[] createIncrementVector(double m, double n, double increment){
		
		if (m>n){
			throw new ArithmeticException(); // for now, can extend later if needed
		}		
		
		double[] out = new double[(int) Math.floor((n-m)/increment) + 1];
		double k = 0;
		for (int i=0;i<out.length;i++){
			out[i] = k+m;
			k+=increment;
		}
		return out;
	}
	
	// Interpolate Matrix
	public static double[][] interpolateMatrix(double[][] matrix, int newsize){
		// preallocate new matrix
		double[][] matrixOut = new double[newsize][matrix[0].length];
		// assign columns
		for (int i = 0; i < matrixOut[0].length; i++) {
		    matrixOut = MyUtilities.assignColumn(matrixOut, MyUtilities.interpolateVector(MyUtilities.fetchColumn(matrix, i), newsize), i);
		}
		return matrixOut;
	}
	
	// Fetch column from Matrix
	public static double[] fetchColumn(double[][] matrix, int columnIndex){
		double[] out = new double[matrix.length];
		for (int i=0;i<matrix.length;i++){
			out[i] = matrix[i][columnIndex];
		}
		return out;
	}
	
	// compare double arrays
	public static boolean compareDoubleArrays(double[] a, double[] b){
		if (a.length != b.length){
			return false;
		}
		double tol = 0.000001;
		for (int k=0;k<a.length;k++){
			if (Math.abs(a[k]-b[k])>tol){
				return false;
			}
		}
		return true;
	}
	
	// compare double matrices
	public static boolean compareDoubleMatrices(double[][] a, double[][] b){
		if (a.length != b.length){
			return false;
		}
		double tol = 0.000001;
		for (int i=0;i<a.length;i++){
			if (!MyUtilities.compareDoubleArrays(a[i], b[i])){
				return false;
			}
		}
		return true;
	}
	
	// remove column from matrix
	public static double[][] removeColumn(double[][] matrix, int columnIndex){
		double[][] matrixOut = new double[matrix.length][matrix[0].length-1];
		int k = 0;
		columnIndex--;
		while (k<matrixOut[0].length){
			if (k<columnIndex){
				matrixOut = MyUtilities.assignColumn(matrixOut, MyUtilities.fetchColumn(matrix, k), k);
			} else {
				matrixOut = MyUtilities.assignColumn(matrixOut, MyUtilities.fetchColumn(matrix, k+1), k);
			}
			k++;
		}
		return matrixOut;
	}
	
	// matrix addition
	public static double[][] addMatrices(double[][] M1, double[][] M2){
		if (M1.length != M2.length | M1[0].length != M2[0].length){
			throw new ArithmeticException();
		}
		double[][] out = new double[M1.length][M1[0].length];
		for (int i = 0;i<M1.length;i++){
			out[i] = addVectors(M1[i], M2[i]);
		}
		return out;
	}
	
	// matrix absolute value
	public static double[][] matrixAbsValue(double[][] M){
		double[][] out = new double[M.length][M[0].length];
		for (int i=0;i<M.length;i++){
			for (int j=0;j<M[i].length;j++){
				out[i][j] = Math.abs(M[i][j]);
			}
			
		}
		return out;
	}	
	
	
	
	// TODO: Write tests for the methods below this
	// create vector of ones
	public static double[] onesVector(int size){
		double[] out = new double[size];
		Arrays.fill(out, 1);
		return out;
	}
	
	// create matrix of zeros
	public static double[][] zerosMatrix(int size1,int size2){
		double[][] out = new double[size1][size2];
		for(int i = 0; i < out.length; i++) {
		    Arrays.fill(out[i], 0);
		}
		return out;
	}
	
	// assign vector to a given row of a matrix
	public static double[][] assignRow(double[][] matrix, double[] v, int index){
		for (int i=0;i<matrix[index].length;i++){
			matrix[index][i] = v[i];
		}
		return matrix;
	}
	
    // element-wise subtraction of two vectors
	public static double[] subtractVectors(double[] v1, double[] v2) {
		if (v1.length != v2.length){
			throw new ArithmeticException();
		}
		double[] out = new double[v1.length];
		
		for(int i = 0;i<v1.length;i++){
			
			out[i] = v1[i] - v2[i];
			
		}
		
		return out;
	}
	
	// vector mean
	public static double meanVector(double[] v){
		double sum = 0;
		for (int i=0;i<v.length;i++){
			sum += v[i];
		}
		return sum / v.length;
	}
	
	// matrix columnwise means
	public static double[] meanColumns(double[][] M){
		double[] out = new double[M[0].length];
		for (int i=0;i<M[0].length;i++){
			out[i] = meanVector(fetchColumn(M,i+1));
		}
		return out;
	}	
	
	// Fetch row from Matrix
	public static double[] fetchRow(double[][] matrix, int rowIndex){
		double[] out = new double[matrix.length];
		for (int i=0;i<matrix[0].length;i++){
			out[i] = matrix[rowIndex-1][i];
		}
		return out;
	}

}
