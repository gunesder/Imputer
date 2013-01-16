package edu.berkeley.path.imputer;

import org.apache.commons.math.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math.analysis.polynomials.PolynomialSplineFunction;
import java.lang.*;

public class MyUtilities {
	
	// *********************************************************************************************
	// ************ Array operations ***************************************************************
	// *********************************************************************************************
	
	// Scaling
	public static double[] scaleVector(double[] d, double s){
		for (int k=0;k<d.length;k++){
			d[k]*=s;			
		}
		return d;
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
			matrix[k][columnIndex-1] = vector[k];			
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
		    matrixOut = MyUtilities.assignColumn(matrixOut, MyUtilities.interpolateVector(MyUtilities.fetchColumn(matrix, i+1), newsize), i+1);
		}
		return matrixOut;
	}
	
	// Fetch column from Matrix
	public static double[] fetchColumn(double[][] matrix, int columnIndex){
		double[] out = new double[matrix.length];
		for (int i=0;i<matrix.length;i++){
			out[i] = matrix[i][columnIndex-1];
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
				matrixOut = MyUtilities.assignColumn(matrixOut, MyUtilities.fetchColumn(matrix, k+1), k+1);
			} else {
				matrixOut = MyUtilities.assignColumn(matrixOut, MyUtilities.fetchColumn(matrix, k+2), k+1);
			}
			k++;
		}
		return matrixOut;
	}

}
