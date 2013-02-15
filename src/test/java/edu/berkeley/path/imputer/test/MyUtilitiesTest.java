package edu.berkeley.path.imputer.test;

import java.util.Arrays;

import edu.berkeley.path.imputer.MyUtilities;
import junit.framework.TestCase;

public class MyUtilitiesTest extends TestCase {

	public void testScaleVector() {
		double[] expected = {0.5,1,1.5};
		double[] input = {1,2,3};
		assertTrue(Arrays.equals(expected, MyUtilities.scaleVector(input, 0.5)));
	}
	
	public void testScaleMatrix(){
		double[][] input = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
		double[][] expected = {{2,4,6,8,10},{12,14,16,18,20},{22,24,26,28,30}};
		assertTrue(MyUtilities.compareDoubleMatrices(expected, MyUtilities.scaleMatrix(input, 2)));
	}

	public void testAddVectors() {
		double[] expected = {0,12,3.5};
		double[] input1 = {1,2,3};
		double[] input2 = {-1,10,0.5};
		assertTrue(Arrays.equals(expected, MyUtilities.addVectors(input1, input2)));
	}

	public void testAssignColumn() {
		// fat matrix
		double[][] expected = {{1,31,3,4,5},{6,31,8,9,10},{11,31,13,14,15}};
		double[][] input1 = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
		double[] input2 = {31,31,31};
		assertTrue(Arrays.deepEquals(expected, MyUtilities.assignColumn(input1, input2, 1)));
		// tall matrix
		double[][] expected2 = {{1,31,3},{4,31,6},{7,31,9},{10,31,12},{13,31,15}};
		double[][] input3 = {{1,2,3},{4,5,6},{7,8,9},{10,11,12},{13,14,15}};
		double[] input4 = {31,31,31,31,31};
		assertTrue(Arrays.deepEquals(expected2, MyUtilities.assignColumn(input3, input4, 1)));		
	}

	public void testAppendColumn() {
		// fat matrix
		double[][] expected = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
		double[][] input1 = {{1,2,3,4},{6,7,8,9},{11,12,13,14}};
		double[] input2 = {5,10,15};
		assertTrue(Arrays.deepEquals(expected, MyUtilities.appendColumn(input1, input2)));
		// tall matrix
		double[][] expected2 = {{1,2,3},{4,5,6},{7,8,9},{10,11,12},{13,14,15}};
		double[][] input3 = {{1,2},{4,5},{7,8},{10,11},{13,14}};
		double[] input4 = {3,6,9,12,15};
		assertTrue(Arrays.deepEquals(expected2, MyUtilities.appendColumn(input3, input4)));		
	}

	public void testSumColumns() {
		// fat matrix
		double[] expected = {10,30,50};
		double[][] input1 = {{1,2,3,4},{6,7,8,9},{11,12,13,14}};
		assertTrue(Arrays.equals(expected, MyUtilities.sumColumns(input1)));
		// tall matrix
		double[] expected2 = {3,9,15,21,27};
		double[][] input2 = {{1,2},{4,5},{7,8},{10,11},{13,14}};
		assertTrue(Arrays.equals(expected2, MyUtilities.sumColumns(input2))); 
	}

	public void testInterpolateVector() {
		double[] expected = {20,30,40,50,60,70,80};
		double[] input = {20,40,60,80};
		assertTrue(MyUtilities.compareDoubleArrays(expected, MyUtilities.interpolateVector(input,7)));
	}
	
	public void testCreateIncrementVector() {
		double[] expected = {-1,-0.7,-0.4,-0.1,0.2,0.5};
		double m = -1, n = 0.5;
		assertTrue(MyUtilities.compareDoubleArrays(expected, MyUtilities.createIncrementVector(m, n, 0.3)));
	}
	
	public void testFetchColumn(){
		// fat matrix
		double[][] input = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
		double[] expected = {4,9,14};
		assertTrue(Arrays.equals(expected, MyUtilities.fetchColumn(input, 3)));
		// tall matrix
		double[][] input2 = {{1,2,3},{4,5,6},{7,8,9},{10,11,12},{13,14,15}};
		double[] expected2 = {3,6,9,12,15};
		assertTrue(Arrays.equals(expected2, MyUtilities.fetchColumn(input2, 2)));
	}
	
	public void testCompareDoubleMatrices(){
		double[][] input1 = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
		double[][] input2 = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
		assertTrue(MyUtilities.compareDoubleMatrices(input1, input2));
	}
	
	public void testInterpolateMatrix(){
		double[][] input = {{20,120,1},{40,130,2},{60,140,3},{80,150,4}};
		double[][] expected = {{20,120,1},{30,125,1.5},{40,130,2},{50,135,2.5},{60,140,3},{70,145,3.5},{80,150,4}};
		assertTrue(MyUtilities.compareDoubleMatrices(expected,MyUtilities.interpolateMatrix(input, 7)));
	}
	
	public void testRemoveColumn(){
		double[][] input = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
		double[][] expected1 = {{2,3,4,5},{7,8,9,10},{12,13,14,15}};
		double[][] expected2 = {{1,2,3,4},{6,7,8,9},{11,12,13,14}};
		double[][] expected3 = {{1,2,4,5},{6,7,9,10},{11,12,14,15}};
		assertTrue(MyUtilities.compareDoubleMatrices(expected1,MyUtilities.removeColumn(input, 1)));
		assertTrue(MyUtilities.compareDoubleMatrices(expected2,MyUtilities.removeColumn(input, 5)));
		assertTrue(MyUtilities.compareDoubleMatrices(expected3,MyUtilities.removeColumn(input, 3)));
	}
	
	public void testAddMatrices(){
		double[][] input1 = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
		double[][] input2 = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
		double[][] expected = {{2,4,6,8,10},{12,14,16,18,20},{22,24,26,28,30}};
		assertTrue(MyUtilities.compareDoubleMatrices(expected, MyUtilities.addMatrices(input1, input2)));
	}
	
	public void testMatrixAbsValue(){
		double[][] input = {{-1,2,3,-4,5},{6,7,-8,9,10},{11,-12,-13,-14,-15}};
		double[][] expected = {{1,2,3,4,5},{6,7,8,9,10},{11,12,13,14,15}};
		assertTrue(MyUtilities.compareDoubleMatrices(expected, MyUtilities.matrixAbsValue(input)));		
	}
	
	public void testFilter(){
		double[] x = MyUtilities.createIncrementVector(1, 4, 0.2);
		double[] b = MyUtilities.scaleVector(MyUtilities.onesVector(5), 0.2);
		double[] a = {1};
		double[] expected = {0.2,0.44,0.72,1.04,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6};
		assertTrue(MyUtilities.compareDoubleArrays(expected, MyUtilities.filter(b, a, x)));
	}

}
