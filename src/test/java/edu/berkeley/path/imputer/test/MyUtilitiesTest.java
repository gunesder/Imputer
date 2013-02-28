package edu.berkeley.path.imputer.test;

import java.util.Arrays;

import edu.berkeley.path.imputer.MyUtilities;
import junit.framework.TestCase;

public class MyUtilitiesTest extends TestCase {

	public void testScaleVector() {
		Double[] expected = {0.5,1.0,1.5};
		Double[] input = {1.0,2.0,3.0};
		assertTrue(Arrays.equals(expected, MyUtilities.scaleVector(input, 0.5)));
	}
	
	public void testScaleMatrix(){
		Double[][] input = {{1.0,2.0,3.0,4.0,5.0},{6.0,7.0,8.0,9.0,10.0},{11.0,12.0,13.0,14.0,15.0}};
		Double[][] expected = {{2.0,4.0,6.0,8.0,10.0},{12.0,14.0,16.0,18.0,20.0},{22.0,24.0,26.0,28.0,30.0}};
		assertTrue(MyUtilities.compareDoubleMatrices(expected, MyUtilities.scaleMatrix(input, 2.0)));
	}

	public void testAddVectors() {
		Double[] expected = {0.0,12.0,3.5};
		Double[] input1 = {1.0,2.0,3.0};
		Double[] input2 = {-1.0,10.0,0.5};
		assertTrue(Arrays.equals(expected, MyUtilities.addVectors(input1, input2)));
	}

	public void testAssignColumn() {
		// fat matrix
		Double[][] expected = {{1.0,31.0,3.0,4.0,5.0},{6.0,31.0,8.0,9.0,10.0},{11.0,31.0,13.0,14.0,15.0}};
		Double[][] input1 = {{1.0,2.0,3.0,4.0,5.0},{6.0,7.0,8.0,9.0,10.0},{11.0,12.0,13.0,14.0,15.0}};
		Double[] input2 = {31.0,31.0,31.0};
		assertTrue(Arrays.deepEquals(expected, MyUtilities.assignColumn(input1, input2, 1)));
		// tall matrix
		Double[][] expected2 = {{1.0,31.0,3.0},{4.0,31.0,6.0},{7.0,31.0,9.0},{10.0,31.0,12.0},{13.0,31.0,15.0}};
		Double[][] input3 = {{1.0,2.0,3.0},{4.0,5.0,6.0},{7.0,8.0,9.0},{10.0,11.0,12.0},{13.0,14.0,15.0}};
		Double[] input4 = {31.0,31.0,31.0,31.0,31.0};
		assertTrue(Arrays.deepEquals(expected2, MyUtilities.assignColumn(input3, input4, 1)));		
	}

	public void testAppendColumn() {
		// fat matrix
		Double[][] expected = {{1.0,2.0,3.0,4.0,5.0},{6.0,7.0,8.0,9.0,10.0},{11.0,12.0,13.0,14.0,15.0}};
		Double[][] input1 = {{1.0,2.0,3.0,4.0},{6.0,7.0,8.0,9.0},{11.0,12.0,13.0,14.0}};
		Double[] input2 = {5.0,10.0,15.0};
		assertTrue(Arrays.deepEquals(expected, MyUtilities.appendColumn(input1, input2)));
		// tall matrix
		Double[][] expected2 = {{1.0,2.0,3.0},{4.0,5.0,6.0},{7.0,8.0,9.0},{10.0,11.0,12.0},{13.0,14.0,15.0}};
		Double[][] input3 = {{1.0,2.0},{4.0,5.0},{7.0,8.0},{10.0,11.0},{13.0,14.0}};
		Double[] input4 = {3.0,6.0,9.0,12.0,15.0};
		assertTrue(Arrays.deepEquals(expected2, MyUtilities.appendColumn(input3, input4)));		
	}

	public void testSumColumns() {
		// fat matrix
		Double[] expected = {10.0,30.0,50.0};
		Double[][] input1 = {{1.0,2.0,3.0,4.0},{6.0,7.0,8.0,9.0},{11.0,12.0,13.0,14.0}};
		assertTrue(MyUtilities.compareDoubleArrays(expected, MyUtilities.sumColumns(input1)));
		// tall matrix
		Double[] expected2 = {3.0,9.0,15.0,21.0,27.0};
		Double[][] input2 = {{1.0,2.0},{4.0,5.0},{7.0,8.0},{10.0,11.0},{13.0,14.0}};
		assertTrue(MyUtilities.compareDoubleArrays(expected2, MyUtilities.sumColumns(input2))); 
	}

	public void testInterpolateVector() {
		Double[] expected = {20.0,30.0,40.0,50.0,60.0,70.0,80.0};
		Double[] input = {20.0,40.0,60.0,80.0};
		assertTrue(MyUtilities.compareDoubleArrays(expected, MyUtilities.interpolateVector(input,7)));
	}
	
	public void testCreateIncrementVector() {
		Double[] expected = {-1.0,-0.7,-0.4,-0.1,0.2,0.5};
		Double m = -1.0, n = 0.5;
		assertTrue(MyUtilities.compareDoubleArrays(expected, MyUtilities.createIncrementVector(m, n, 0.3)));
	}
	
	public void testFetchColumn(){
		// fat matrix
		Double[][] input = {{1.0,2.0,3.0,4.0,5.0},{6.0,7.0,8.0,9.0,10.0},{11.0,12.0,13.0,14.0,15.0}};
		Double[] expected = {4.0,9.0,14.0};
		assertTrue(Arrays.equals(expected, MyUtilities.fetchColumn(input, 3)));
		// tall matrix
		Double[][] input2 = {{1.0,2.0,3.0},{4.0,5.0,6.0},{7.0,8.0,9.0},{10.0,11.0,12.0},{13.0,14.0,15.0}};
		Double[] expected2 = {3.0,6.0,9.0,12.0,15.0};
		assertTrue(Arrays.equals(expected2, MyUtilities.fetchColumn(input2, 2)));
	}
	
	public void testCompareDoubleMatrices(){
		Double[][] input1 = {{1.0,2.0,3.0,4.0,5.0},{6.0,7.0,8.0,9.0,10.0},{11.0,12.0,13.0,14.0,15.0}};
		Double[][] input2 = {{1.0,2.0,3.0,4.0,5.0},{6.0,7.0,8.0,9.0,10.0},{11.0,12.0,13.0,14.0,15.0}};
		assertTrue(MyUtilities.compareDoubleMatrices(input1, input2));
	}
	
	public void testInterpolateMatrix(){
		Double[][] input = {{20.0,120.0,1.0},{40.0,130.0,2.0},{60.0,140.0,3.0},{80.0,150.0,4.0}};
		Double[][] expected = {{20.0,120.0,1.0},{30.0,125.0,1.5},{40.0,130.0,2.0},{50.0,135.0,2.5},{60.0,140.0,3.0},{70.0,145.0,3.5},{80.0,150.0,4.0}};
		assertTrue(MyUtilities.compareDoubleMatrices(expected,MyUtilities.interpolateMatrix(input, 7)));
	}
	
	public void testRemoveColumn(){
		Double[][] input = {{1.0,2.0,3.0,4.0,5.0},{6.0,7.0,8.0,9.0,10.0},{11.0,12.0,13.0,14.0,15.0}};
		Double[][] expected1 = {{2.0,3.0,4.0,5.0},{7.0,8.0,9.0,10.0},{12.0,13.0,14.0,15.0}};
		Double[][] expected2 = {{1.0,2.0,3.0,4.0},{6.0,7.0,8.0,9.0},{11.0,12.0,13.0,14.0}};
		Double[][] expected3 = {{1.0,2.0,4.0,5.0},{6.0,7.0,9.0,10.0},{11.0,12.0,14.0,15.0}};
		assertTrue(MyUtilities.compareDoubleMatrices(expected1,MyUtilities.removeColumn(input, 1)));
		assertTrue(MyUtilities.compareDoubleMatrices(expected2,MyUtilities.removeColumn(input, 5)));
		assertTrue(MyUtilities.compareDoubleMatrices(expected3,MyUtilities.removeColumn(input, 3)));
	}
	
	public void testAddMatrices(){
		Double[][] input1 = {{1.0,2.0,3.0,4.0,5.0},{6.0,7.0,8.0,9.0,10.0},{11.0,12.0,13.0,14.0,15.0}};
		Double[][] input2 = {{1.0,2.0,3.0,4.0,5.0},{6.0,7.0,8.0,9.0,10.0},{11.0,12.0,13.0,14.0,15.0}};
		Double[][] expected = {{2.0,4.0,6.0,8.0,10.0},{12.0,14.0,16.0,18.0,20.0},{22.0,24.0,26.0,28.0,30.0}};
		assertTrue(MyUtilities.compareDoubleMatrices(expected, MyUtilities.addMatrices(input1, input2)));
	}
	
	public void testMatrixAbsValue(){
		Double[][] input = {{-1.0,2.0,3.0,-4.0,5.0},{6.0,7.0,-8.0,9.0,10.0},{11.0,-12.0,-13.0,-14.0,-15.0}};
		Double[][] expected = {{1.0,2.0,3.0,4.0,5.0},{6.0,7.0,8.0,9.0,10.0},{11.0,12.0,13.0,14.0,15.0}};
		assertTrue(MyUtilities.compareDoubleMatrices(expected, MyUtilities.matrixAbsValue(input)));		
	}
	
	public void testFilter(){
		Double[] x = MyUtilities.createIncrementVector(1.0, 4.0, 0.2);
		Double[] b = MyUtilities.scaleVector(MyUtilities.onesVector(5), 0.2);
		Double[] a = {1.0};
		Double[] expected = {0.2,0.44,0.72,1.04,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6};
		assertTrue(MyUtilities.compareDoubleArrays(expected, MyUtilities.filter(b, a, x)));
	}
	
	public void testMeanVector(){
		Double[] a = {-2.0,3.0,6.0,-18.0,72.0};
		assertTrue(MyUtilities.meanVector(a) > 12.2-0.000001 && MyUtilities.meanVector(a) < 12.2+0.000001);
	}
	
	public void testMeanColumns(){
		Double[][] input1 = {{1.0,2.0,3.0,4.0,5.0},{6.0,7.0,8.0,9.0,10.0},{11.0,12.0,13.0,14.0,15.0}};
		Double[] expected = {6.0,7.0,8.0,9.0,10.0};
		assertTrue(MyUtilities.compareDoubleArrays(expected, MyUtilities.meanColumns(input1)));
	}
	
	public void testMeanRows(){
		Double[][] input1 = {{1.0,2.0,3.0,4.0,5.0},{6.0,7.0,8.0,9.0,10.0},{11.0,12.0,13.0,14.0,15.0}};
		Double[] expected = {3.0,8.0,13.0};
		assertTrue(MyUtilities.compareDoubleArrays(expected, MyUtilities.meanRows(input1)));
	}
	
	public void testMultiplyMatrices(){
		Double[][] input1 = {{1.0,2.0,3.0},{4.0,5.0,6.0},{7.0,8.0,9.0}};
		Double[][] input2 = {{1.0,2.0},{3.0,4.0},{5.0,6.0}};
		Double[][] expected = {{22.0,28.0},{49.0,64.0},{76.0,100.0}};
		assertTrue(MyUtilities.compareDoubleMatrices(expected, MyUtilities.multiplyMatrices(input1,input2)));
	}

}
