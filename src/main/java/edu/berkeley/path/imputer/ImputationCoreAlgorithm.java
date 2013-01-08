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
		private ArrayList<Boolean> imputeOR = new ArrayList<Boolean>();
		private ArrayList<Boolean> imputeFR = new ArrayList<Boolean>();
		private ArrayList<Boolean> orPresent = new ArrayList<Boolean>();
		private ArrayList<Boolean> frPresent = new ArrayList<Boolean>();
		// derived fields (matrix quantities)
		private BlockRealMatrix measuredDensity = new BlockRealMatrix(detectorList.get(null).getFlowData().size(),cellData.size());
		private BlockRealMatrix measuredFlow = new BlockRealMatrix(detectorList.get(null).getFlowData().size(),cellData.size());
		private BlockRealMatrix measuredSpeed = new BlockRealMatrix(detectorList.get(null).getFlowData().size(),cellData.size());
		// learning algorithm parameters
		// TODO: learning algo parameters
		
		
	
	// getters and setters
		// TODO: generate getters and setters
	
	// constructors
	public ImputationCoreAlgorithm(LinkedList<Cell> cells, HashMap<Integer,Detector> detectors){
		this.cellData = cells;
		this.detectorList = detectors;
	}
	
	public ImputationCoreAlgorithm(LinkedList<Cell> cells, HashMap<Integer,Detector> detectors, double simTimeStep){
		this.cellData = cells;
		this.detectorList = detectors;
		this.simulationTimeStep = simTimeStep;
	}
	
	public ImputationCoreAlgorithm(LinkedList<Cell> cells, HashMap<Integer,Detector> detectors, double simTimeStep, double demTimeStep){
		this.cellData = cells;
		this.detectorList = detectors;
		this.simulationTimeStep = simTimeStep;
		this.demandTimeStep = demTimeStep;
	}
	
	// methods
	public void run(){
		
		// double firstPM = cellData.getFirst().getDetectorML().getAbsolutePM();
		// fill parameter vectors, data matrices and interpolate to appropriate simulation time step
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
			
			
			
		}
	}

}
