package edu.berkeley.path.imputer;

import org.apache.commons.math.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math.analysis.polynomials.PolynomialSplineFunction;

import scala.actors.threadpool.Arrays;

import java.io.File;
import java.io.IOException;
import java.lang.*;
import java.util.ArrayList;

import jxl.NumberCell;
import jxl.Workbook;
import jxl.read.biff.BiffException;

public class MyUtilities {
	
	// *********************************************************************************************
	// ************ Array operations ***************************************************************
	// *********************************************************************************************
	
	// Scale Vector
	public static double[] scaleVector(double[] d, double s){
		double[] out = new double[d.length];
		for (int k=0;k<d.length;k++){
			out[k] = d[k]*s;			
		}
		return out;
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
		double[] evaluationPoints = MyUtilities.createIncrementVector(0,(double) vector.length-1,(double) (vector.length-1)/(newsize-1));
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
			out[i] = meanVector(fetchColumn(M,i));
		}
		return out;
	}	
	
	// matrix rowwise means
	public static double[] meanRows(double[][] M){
		double[] out = new double[M.length];
		for (int i=0;i<M.length;i++){
			out[i] = meanVector(fetchRow(M,i));
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
	
	// Fetch row from Matrix
	public static double[] fetchRow(double[][] matrix, int rowIndex){
		double[] out = new double[matrix[0].length];
		for (int i=0;i<matrix[0].length;i++){
			out[i] = matrix[rowIndex][i];
		}
		return out;
	}
	
	// Savitzky - Golay smoothing (can overload with other smoothing algorithms if required)
	public static double[] smoothFWM(double[] data, int f, int k){
		
		/* (x,y) are given data. f is the frame length to be taken, should
		be an odd number. k is the degree of polynomial filter. It should
		be less than f.

		Reference: Orfanidis, S.J., Introduction to Signal Processing,
		Prentice-Hall, Englewood Cliffs, NJ, 1996. */
		
		double[] x = MyUtilities.createIncrementVector(1, data.length, 1);
		int n = x.length;
		f = (int) Math.floor(f);
		f = Math.min(f, n);
		f = f - (f-1)%2; // will subtract 1 if frame is even.
		double[] diffx = MyUtilities.onesVector(n-1);
		// skipping checks for NaNs and Nulls in the data. Also skipping check if span is less than degree
		
		int hf = (f-1)/2; // half frame length
		
		double[][] v = MyUtilities.zerosMatrix(f, k+1);
		for (int i=0;i<v.length;i++){
			for (int j=0;j<v[0].length;j++){
				v[i][j]++;
			}
		}
		
		double[] t = MyUtilities.createIncrementVector(-hf, hf, 1);
		double[] dummy = new double[t.length];
		for (int i=0;i<k+1;i++){
			
			dummy = MyUtilities.onesVector(dummy.length);
			for (int j=0;j<dummy.length;j++){
				dummy[j] = Math.pow(t[j],i);
			}
		
			MyUtilities.assignColumn(v, dummy, i);
		}
		
		QRDecomposer decomposer = new QRDecomposer(v);
		double[][] q = decomposer.getQ();
		
		double[][] dummy1 = new double[2][1]; // this is my ad hoc way of defining vectors as 2d arrays. This is a bit of a hack to conform with the matrix multiplication code I took off the web. 
		dummy1[0] = MyUtilities.fetchRow(q, hf);
		dummy1 = MyUtilities.multiplyMatrices(q, dummy1);
		
		// ymid
		double[] b = MyUtilities.fetchColumn(dummy1, 0);
		double[] a = {1};
		double[] ymid = MyUtilities.filter(b, a, data);
		
		// ybegin
		double[][] dummy2 = new double[hf][q[0].length];
		for (int i=0;i<dummy2.length;i++){
			dummy2[i] = MyUtilities.fetchRow(q, i);
		}
		
		double[][] dummy3 = MyUtilities.multiplyMatrices(dummy2, MyUtilities.transposeMatrix(q));
		double[][] ybegin = new double[2][1];
		ybegin[0] = new double[f];
		for (int i=0;i<ybegin[0].length;i++){
			ybegin[0][i] = data[i];
		}
		ybegin = MyUtilities.multiplyMatrices(dummy3, ybegin);
		
		// yend
		double[][] dummy4 = new double[hf][q[0].length];
		for (int i=hf+1;i<q.length;i++){
			dummy4 = MyUtilities.assignRow(dummy4, MyUtilities.fetchRow(q, i), i-hf-1);
		}
		double[][] dummy5 = MyUtilities.multiplyMatrices(dummy4, MyUtilities.transposeMatrix(q));
		double[][] yend = new double[2][1];
		yend[0] = new double[f];
		for (int i=data.length-f;i<data.length;i++){
			yend[0][i-(data.length-f)] = data[i];
		}
		yend = MyUtilities.multiplyMatrices(dummy5, yend);
		
		// putting the y's together for the output
		double[] out = new double[data.length];
		for (int i=0;i<out.length;i++){
			if (i<ybegin.length){
				out[i] = ybegin[i][0];
			} else if (i>=out.length-yend.length){
				out[i] = yend[i-(out.length-yend.length)][0];
			} else {
				out[i] = ymid[i];
			}
		}
		return out;		
	
		// non-uniform x is ignored because it's never necessary for the cases needed for imputation
		
	}
	
	// implementation of the Matlab filter function
	public static double[] filter(double[] b, double[] a, double[] data){
		
		/* b is the vector of coefficients of the feedforward filter and a is the 
		 * vector of coefficients for the feedback filter. Refer to Matlab documentation 
		 * for the function "filter" for further detail.
		 * 
		 * Reference:  Oppenheim, A. V. and R.W. Schafer. Discrete-Time Signal Processing,
		 * Englewood Cliffs, NJ: Prentice-Hall, 1989, pp. 311-312. */
		
		double[] out = new double[data.length];
		
		for (int i=0;i<data.length;i++){
			// discrete time filtering, difference equation consists of feedforward terms - feedback terms.
			
			// positive term (feedforward):
			double posTerm = 0;
			int k = i;
			for (int j=0;j<b.length;j++){
				if (k>=0){
					posTerm += b[j]*data[k];
					k--;
				} else {
					break;
				}
			}
			
			// negative term (feedback):
			double negTerm = 0;
			if (a.length>1){
				k = i-1;
				for (int j=1;j<a.length;j++){
					if (k>=0){
						negTerm += a[j]*data[k];
						k--;
					} else {
						break;
					}
				}
			}
			
			out[i] = posTerm - negTerm;
			
			
		}
		
		return out;
		
	}
	
	// Matrix multiplication
	public static double[][] multiplyMatrices(double a[][], double b[][]) {
		   
		int aRows = a.length, aColumns = a[0].length, bRows = b.length, bColumns = b[0].length;
		
		// modify dimensions for vector inputs disguised as 2d arrays
		// a is a row vector:
		if (a[0].length>a[1].length){
			aColumns = Math.max(aRows, aColumns);
			aRows = 1;
		}
		// b is a column vector
		if (b[0].length>b[1].length){
			bRows = Math.max(bRows, bColumns);
			bColumns = 1;
			double[][] dummy = new double[b[0].length][1]; 
			for (int i=0;i<b[0].length;i++){
				dummy[i][0] = b[0][i];
			}
			b = dummy;
		}
		
		   
		if ( aColumns != bRows ) {
			throw new IllegalArgumentException("A:Rows: " + aColumns + " did not match B:Columns " + bRows + ".");
		}
		   
		double[][] result = new double[aRows][bColumns];
		   
		for(int i = 0; i < aRows; i++) { // aRow
			for(int j = 0; j < bColumns; j++) { // bColumn
				for(int k = 0; k < aColumns; k++) { // aColumn
					result[i][j] += a[i][k] * b[k][j];
		        }
		    }  
		}
		   
		return result;
	}
	
	// Matrix transpose
	public static double[][] transposeMatrix(double[][] A){
		
		double[][] out = new double[A[0].length][A.length];
		
		for (int i = 0;i<out.length;i++){
			out[i] = MyUtilities.fetchColumn(A, i);
		}
		
		return out;
		
	}
	
	public static double[][] read2DArrayFromExcel(String filename) throws BiffException, IOException{
		
		Workbook workbook = Workbook.getWorkbook(new File(filename));
		double[][] out = new double[workbook.getSheet(0).getRows()][workbook.getSheet(0).getColumns()];
		for (int i=0;i<workbook.getSheet(0).getRows();i++){
			for (int j=0;j<workbook.getSheet(0).getColumns();j++){
				NumberCell nc = (NumberCell) workbook.getSheet(0).getCell(j, i); 
				out[i][j] = (double) nc.getValue();
			}
		}
		
		return out;	
		
	}
	
	public static double[][] reshapeVectorIntoMatrix(double[] v, int m, int n){
		double[][] out = new double[m][n];
		int k=0;
		for (int i=0;i<m;i++){
			for (int j=0;j<n;j++){
				out[i][j] = v[k++];
			}
		}
		return out;
	}
		
	

}
