package edu.berkeley.path.imputer;

import java.util.*;
import org.apache.commons.math.linear.BlockRealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.analysis.interpolation.LinearInterpolator;

public class ImputationCoreAlgorithm {
	
	// fields
		// primary fields
		private LinkedList<Cell> cellData = new LinkedList<Cell>();
		private HashMap<Integer,Detector> detectorList = new HashMap<Integer,Detector>();
		private double demandTimeStep = 5/60; // [hours] 5 minute default
		private double simulationTimeStep = 5/60/60; // [hours] 5 second default
		// derived fields (vector quantities)
		private ArrayList<Double> linkLength = new ArrayList<Double>();
		private ArrayList<Integer> lane = new ArrayList<Integer>();
		private ArrayList<Double> vf = new ArrayList<Double>();
		private ArrayList<Double> w = new ArrayList<Double>();
		private ArrayList<Double> rhojam = new ArrayList<Double>();
		private ArrayList<Double> qmax = new ArrayList<Double>();
		private ArrayList<Double> maxBeta = new ArrayList<Double>();
		private ArrayList<Double> maxD = new ArrayList<Double>();
		private ArrayList<Boolean> imputeOR = new ArrayList<Boolean>();
		private ArrayList<Boolean> imputeFR = new ArrayList<Boolean>();
		private ArrayList<Boolean> orPresent = new ArrayList<Boolean>();
		private ArrayList<Boolean> frPresent = new ArrayList<Boolean>();
		// derived fields (matrix quantities)
		private BlockRealMatrix measuredDensity = new BlockRealMatrix(288,1);
		private BlockRealMatrix measuredFlow = new BlockRealMatrix(288,1);
		private BlockRealMatrix measuredSpeed = new BlockRealMatrix(288,1);
		private BlockRealMatrix measuredOrFlow = new BlockRealMatrix(288,1);
		private BlockRealMatrix measuredFrFlow = new BlockRealMatrix(288,1);
		// learning algorithm parameters
		// TODO: learning algo parameters		
	
	// getters and setters
		// TODO: generate getters and setters
	
	// constructors
	public ImputationCoreAlgorithm(LinkedList<Cell> cells, HashMap<Integer,Detector> detectors){
		this.cellData = cells;
		this.detectorList = detectors;
		this.initializeDataMatrices();
	}
	
	public ImputationCoreAlgorithm(LinkedList<Cell> cells, HashMap<Integer,Detector> detectors, double simTimeStep){
		this.cellData = cells;
		this.detectorList = detectors;
		this.simulationTimeStep = simTimeStep;
		this.initializeDataMatrices();
	}
	
	public ImputationCoreAlgorithm(LinkedList<Cell> cells, HashMap<Integer,Detector> detectors, double simTimeStep, double demTimeStep){
		this.cellData = cells;
		this.detectorList = detectors;
		this.simulationTimeStep = simTimeStep;
		this.demandTimeStep = demTimeStep;
		this.initializeDataMatrices();
	}
	
	// methods
	public void run(){
		
		// *********************************************************************************************
		// fill parameter vectors, data matrices and interpolate to appropriate simulation time step ***
		// *********************************************************************************************
		
		int i = 0;
		for (Cell c: cellData){
			
			// Lengths
			double dummyLength = 0;
			for(Link l:c.getLinks()){
				dummyLength = dummyLength + l.getLength();
			}
			linkLength.add(dummyLength);
			
			// Lanes
			lane.add(c.getLinks().getFirst().getLanesML()+c.getLinks().getFirst().getLanesHOV());
			
			// Speed Measurements
			RealVector tempSPD_ML = (RealVector) c.getDetectorML().getSpeedData();
			RealVector tempSPD_HOV = (RealVector) c.getDetectorHOV().getSpeedData();
			tempSPD_ML.mapMultiplyToSelf(c.getLinks().getFirst().getLanesML());
			tempSPD_HOV.mapMultiplyToSelf(c.getLinks().getFirst().getLanesHOV());
			RealVector tempSPD = tempSPD_ML.add(tempSPD_HOV);
			measuredSpeed.setColumnVector(i, tempSPD.mapDivideToSelf(c.getLinks().getFirst().getLanesML()+c.getLinks().getFirst().getLanesHOV()));
			
			// Density Measurements
			RealVector tempDTY_ML = (RealVector) c.getDetectorML().getDensityData();
			RealVector tempDTY_HOV = (RealVector) c.getDetectorHOV().getDensityData();
			RealVector tempDTY = tempDTY_ML.add(tempDTY_HOV);
			measuredDensity.setColumnVector(i, tempDTY.mapMultiplyToSelf(dummyLength));
			
			// Flow Measurements
			RealVector tempFLW_ML = (RealVector) c.getDetectorML().getFlowData();
			RealVector tempFLW_HOV = (RealVector) c.getDetectorHOV().getFlowData();
			RealVector tempFLW = tempFLW_ML.add(tempFLW_HOV);
			measuredFlow.setColumnVector(i, tempFLW.mapMultiplyToSelf(this.simulationTimeStep));
			
			// Fundamental Diagram Parameters
			if (c.getDetectorHOV().getSensorID() != 0) { // case where HOV lane detection is separate
				
				// Mainline
				double rhocrit_ML = c.getDetectorML().getFdParams().getRho_crit()*dummyLength;
				double rhojam_ML = c.getDetectorML().getFdParams().getRho_jam()*dummyLength;
				double qmax_ML = c.getDetectorML().getFdParams().getQ_max()*this.simulationTimeStep;
				
				// HOV
				double rhocrit_HOV = c.getDetectorHOV().getFdParams().getRho_crit()*dummyLength;
				double rhojam_HOV = c.getDetectorHOV().getFdParams().getRho_jam()*dummyLength;
				double qmax_HOV = c.getDetectorHOV().getFdParams().getQ_max()*this.simulationTimeStep;
				
				// Aggregate
				vf.add((qmax_ML+qmax_HOV)/(rhocrit_ML+rhocrit_HOV));
				w.add((qmax_ML+qmax_HOV)/(rhojam_ML+rhojam_HOV-rhocrit_ML-rhocrit_HOV));
				rhojam.add(rhojam_ML+rhojam_HOV);
				qmax.add(qmax_ML+qmax_HOV);
			
			} else {
				
				vf.add(c.getDetectorML().getFdParams().getVf()*this.simulationTimeStep/dummyLength);
				w.add(c.getDetectorML().getFdParams().getW()*this.simulationTimeStep/dummyLength);
				rhojam.add(c.getDetectorML().getFdParams().getRho_jam()*dummyLength);
				qmax.add(c.getDetectorML().getFdParams().getQ_max()*this.simulationTimeStep);
				
			}
			
			// Known Ramp Flows
			RealVector dummyVector = c.getMeasuredOnrampFlow().getColumnVector(0);
			for (int j = 1; j<c.getLinks().size();j++){
				dummyVector = dummyVector.add(c.getMeasuredOnrampFlow().getColumnVector(j));
			}
			measuredOrFlow.setColumnVector(i, dummyVector);
			
			RealVector dummyVector2 = c.getMeasuredOfframpFlow().getColumnVector(0);
			for (int j = 1; j<c.getLinks().size();j++){
				dummyVector2 = dummyVector2.add(c.getMeasuredOfframpFlow().getColumnVector(j));
			}
			measuredFrFlow.setColumnVector(i, dummyVector2);
			
			// Ramps present and need imputation boolean vectors
			boolean dummy = false;
			for (boolean b:c.isImputeOR()){
				if (b) dummy = true;
			}
			imputeOR.add(dummy);
			
			boolean dummy2 = false;
			for (boolean b:c.isImputeFR()){
				if (b) dummy2 = true;
			}
			imputeFR.add(dummy2);
			
			boolean dummy3 = false;
			for (int m:c.getOnrampsPerLink()){
				if (m>0) dummy3 = true;
			}
			orPresent.add(dummy3);
			
			int val = dummy3? 1 : 0;
			maxD.add(val*10000*simulationTimeStep); // Hard-coded maximum demand cap
			
			boolean dummy4 = false;
			for (int m:c.getOfframpsPerLink()){
				if (m>0) dummy4 = true;
			}
			frPresent.add(dummy4);			
			
			int val2 = dummy4? 1 : 0;
			maxBeta.add(val2*0.9); // Hard-coded maximum beta
			
		}
		
		measuredSpeed = this.interpolate((int) (24/simulationTimeStep), measuredSpeed);
		measuredFlow = this.interpolate((int) (24/simulationTimeStep), measuredFlow);
		measuredDensity = this.interpolate((int) (24/simulationTimeStep), measuredDensity);
		measuredOrFlow = this.interpolate((int) (24/simulationTimeStep), measuredOrFlow);
		measuredFrFlow = this.interpolate((int) (24/simulationTimeStep), measuredFrFlow);
		
		
		
		
	}
	
	private BlockRealMatrix interpolate(int newsize, BlockRealMatrix matrix){
		// allocate new matrix with new size
		BlockRealMatrix resultMatrix = new BlockRealMatrix(newsize,cellData.size());
		// create array of interpolation points
		double[] interpolationPoints = new double[newsize];
		for (int k=0;k<interpolationPoints.length;k++){
			interpolationPoints[k] = k;
		}
		// interpolate and assing to output
		LinearInterpolator interliPolat = new LinearInterpolator();
		for (int k=0;k<matrix.getColumnDimension();k++){
			resultMatrix.setColumnVector(k, (RealVector) interliPolat.interpolate(interpolationPoints, matrix.getColumn(k)));
		}		
		return resultMatrix;
	}
	
	private void initializeDataMatrices(){
		
		measuredDensity = new BlockRealMatrix(detectorList.values().iterator().next().getDensityData().size(),cellData.size());
		measuredFlow = new BlockRealMatrix(detectorList.values().iterator().next().getDensityData().size(),cellData.size());
		measuredSpeed = new BlockRealMatrix(detectorList.values().iterator().next().getDensityData().size(),cellData.size());
		measuredOrFlow = new BlockRealMatrix(detectorList.values().iterator().next().getDensityData().size(),cellData.size());
		measuredFrFlow = new BlockRealMatrix(detectorList.values().iterator().next().getDensityData().size(),cellData.size());
		
	}

}
