package edu.berkeley.path.imputer;

import java.util.*;

public class ImputationCoreAlgorithm {
	
	// fields
		// primary fields
		private LinkedList<Cell> cellData = new LinkedList<Cell>();
		private HashMap<Integer,Detector> detectorList = new HashMap<Integer,Detector>();
		private double demandTimeStep = 5.0/60.0; // [hours] 5 minute default
		private double simulationTimeStep = 5.0/60.0/60.0; // [hours] 5 second default
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
		private double[][] measuredDensity = new double[289][1];
		private double[][] measuredFlow = new double[289][1];
		private double[][] measuredSpeed = new double[289][1];
		private double[][] measuredOrFlow = new double[289][1];
		private double[][] measuredFrFlow = new double[289][1];
		// learning algorithm parameters
		private boolean downBoundaryCongested = false;
		// TODO: learning algo parameters		
	
	// getters and setters
		// TODO: generate getters and setters
	
	// constructors
	public ImputationCoreAlgorithm(LinkedList<Cell> cells, HashMap<Integer,Detector> detectors){
		this.cellData = cells;
		this.detectorList = detectors;
		this.initializeDataMatrices();
	}
	
	public ImputationCoreAlgorithm(LinkedList<Cell> cells, HashMap<Integer,Detector> detectors, boolean flwBoundary){
		this.cellData = cells;
		this.detectorList = detectors;
		this.downBoundaryCongested = flwBoundary;
		this.initializeDataMatrices();
	}
	
	public ImputationCoreAlgorithm(LinkedList<Cell> cells, HashMap<Integer,Detector> detectors, boolean flwBoundary, double simTimeStep){
		this.cellData = cells;
		this.detectorList = detectors;
		this.downBoundaryCongested = flwBoundary;
		this.simulationTimeStep = simTimeStep;
		this.initializeDataMatrices();
	}
	
	public ImputationCoreAlgorithm(LinkedList<Cell> cells, HashMap<Integer,Detector> detectors, boolean flwBoundary, double simTimeStep, double demTimeStep){
		this.cellData = cells;
		this.detectorList = detectors;
		this.downBoundaryCongested = flwBoundary;
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
			
			// Determine if the cell has separate HOV detection
			boolean separateHOV = false;
			if (c.getDetectorHOV() != null) separateHOV = true;
			// Lengths
			double dummyLength = 0;
			for(Link l:c.getLinks()){
				dummyLength = dummyLength + l.getLength();
			}
			linkLength.add(dummyLength);
			
			// Lanes
			lane.add(c.getLinks().getFirst().getLanesML()+c.getLinks().getFirst().getLanesHOV());
			
			// Speed Measurements
			if (separateHOV){
				double[] tempSPD_ML = c.getDetectorML().getSpeedDataArray();
				double[] tempSPD_HOV = c.getDetectorHOV().getSpeedDataArray();
				tempSPD_ML = MyUtilities.scaleVector(tempSPD_ML,c.getLinks().getFirst().getLanesML());
				tempSPD_HOV = MyUtilities.scaleVector(tempSPD_HOV,c.getLinks().getFirst().getLanesHOV());
				double[] tempSPD = MyUtilities.addVectors(tempSPD_ML, tempSPD_HOV);
				tempSPD = MyUtilities.scaleVector(tempSPD, 1/(c.getLinks().getFirst().getLanesML()+c.getLinks().getFirst().getLanesHOV()));
				measuredSpeed = MyUtilities.assignColumn(measuredSpeed, tempSPD, i+1);
			} else {
				measuredSpeed = MyUtilities.assignColumn(measuredSpeed, c.getDetectorML().getSpeedDataArray(), i+1);
			}
			
			// Density Measurements
			if (separateHOV){
				double[] tempDTY_ML = c.getDetectorML().getDensityDataArray();
				double[] tempDTY_HOV = c.getDetectorHOV().getDensityDataArray();
				double[] tempDTY = MyUtilities.addVectors(tempDTY_ML,tempDTY_HOV);
				measuredDensity = MyUtilities.assignColumn(measuredDensity, MyUtilities.scaleVector(tempDTY,dummyLength*(c.getLinks().getFirst().getLanesML()+c.getLinks().getFirst().getLanesHOV())), i+1);
			} else {
				measuredDensity = MyUtilities.assignColumn(measuredDensity, MyUtilities.scaleVector(c.getDetectorML().getDensityDataArray(),dummyLength*c.getLinks().getFirst().getLanesML()), i+1);
			}
			
			// Flow Measurements
			if (separateHOV){
				double[] tempFLW_ML = c.getDetectorML().getFlowDataArray();
				double[] tempFLW_HOV = c.getDetectorHOV().getFlowDataArray();
				double[] tempFLW = MyUtilities.addVectors(tempFLW_ML,tempFLW_HOV);
				tempFLW = MyUtilities.scaleVector(tempFLW,this.simulationTimeStep*(c.getLinks().getFirst().getLanesML()+c.getLinks().getFirst().getLanesHOV()));
				measuredFlow = MyUtilities.assignColumn(measuredFlow, tempFLW, i+1);
			} else {
				measuredFlow = MyUtilities.assignColumn(measuredFlow, MyUtilities.scaleVector(c.getDetectorML().getFlowDataArray(),this.simulationTimeStep*c.getLinks().getFirst().getLanesML()), i+1);
			}
			
			// Fundamental Diagram Parameters
			if (separateHOV) { // case where HOV lane detection is separate
				
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
			double[] dummyVector = MyUtilities.sumColumns(c.getMeasuredOnrampFlow());
			measuredOrFlow = MyUtilities.assignColumn(measuredOrFlow, dummyVector, i+1);
			
			double[] dummyVector2 = MyUtilities.sumColumns(c.getMeasuredOfframpFlow());
			measuredFrFlow = MyUtilities.assignColumn(measuredFrFlow, dummyVector2, i+1);
			
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
			
			i++;
		}
		
		measuredSpeed = MyUtilities.interpolateMatrix(measuredSpeed, (int) (24/this.simulationTimeStep));
		measuredFlow = MyUtilities.interpolateMatrix(measuredFlow, (int) (24/this.simulationTimeStep));
		measuredDensity = MyUtilities.interpolateMatrix(measuredDensity,(int) (24/this.simulationTimeStep));
		measuredOrFlow = MyUtilities.interpolateMatrix(measuredOrFlow,(int) (24/this.simulationTimeStep));
		measuredFrFlow = MyUtilities.interpolateMatrix(measuredFrFlow,(int) (24/this.simulationTimeStep));
		
		// variables needed for the downstream boundary congested case (line 126 of matlab code)
		double[] boundaryVelocity = MyUtilities.fetchColumn(measuredSpeed, measuredSpeed[0].length);
		double boundaryVF = this.vf.get(this.vf.size()-1)*this.linkLength.get(this.linkLength.size()-1)/this.simulationTimeStep; 
		
		// additional post-processing of variables (line 165 of matlab code)
		double[] orFLW_save = MyUtilities.fetchColumn(measuredOrFlow, 1);
		double[] inputFLW = MyUtilities.addVectors(MyUtilities.fetchColumn(measuredFlow, 1),orFLW_save);
		// skipped the assignment in line 169 of matlab
		
		// n'th cell offramp is paired with the n+1'th cell onramp at the node level. the following manipulations are made to take this into account
		imputeOR.remove(0);
		imputeFR.remove(imputeFR.size()-1);
		orPresent.remove(0);
		frPresent.remove(frPresent.size()-1);
		
		double[][] OrFlow = MyUtilities.removeColumn(measuredOrFlow, 1);
		double[][] FrFlow = MyUtilities.removeColumn(measuredFrFlow, measuredFrFlow[0].length);
		
		double[][] orFlow_Giv = measuredOrFlow; // probably superfluous but keeping up with the matlab code for now
		double[][] frFlow_Giv = measuredFrFlow; // probably superfluous but keeping up with the matlab code for now
		
		ArrayList<Boolean> impute = new ArrayList<Boolean>();
		for (int j=0;i<imputeOR.size();j++){
			impute.add(imputeOR.get(j)|imputeFR.get(j));
		}
		
		int numberOfNodes = measuredDensity[0].length + 1;
		
		double[][] BETA = new double[FrFlow.length][FrFlow[0].length];
		double[][] dj = new double[FrFlow.length][FrFlow[0].length];
		double[] zeroVector = new double[FrFlow.length];
		Arrays.fill(zeroVector, 0);
		BETA = MyUtilities.appendColumn(BETA, zeroVector);
		dj = MyUtilities.appendColumn(dj, zeroVector);
		double[][] OrInp = new double[dj.length][dj[0].length];
		double[] dprev = new double[numberOfNodes-2];
		Arrays.fill(dprev, 0);
		double[][] ErrA = new double[BETA.length][BETA[0].length];
		double[][] ErrB = new double[BETA.length][BETA[0].length];
		
		// skipped line 198 for now
		
		// *********************************************************************************************
		// Derive dj and BETA from known ramp flows ****************************************************
		// *********************************************************************************************
		
		
		
		
		
			
	}
	
	private void initializeDataMatrices(){
		
		measuredDensity = new double[detectorList.values().iterator().next().getDensityData().size()][cellData.size()];
		measuredFlow = new double[detectorList.values().iterator().next().getDensityData().size()][cellData.size()];
		measuredSpeed = new double[detectorList.values().iterator().next().getDensityData().size()][cellData.size()];
		measuredOrFlow = new double[detectorList.values().iterator().next().getDensityData().size()][cellData.size()];
		measuredFrFlow = new double[detectorList.values().iterator().next().getDensityData().size()][cellData.size()];
		
//		Arrays.fill(measuredDensity, 0);
//		Arrays.fill(measuredFlow, 0);
//		Arrays.fill(measuredSpeed, 0);
//		Arrays.fill(measuredOrFlow, 0);
//		Arrays.fill(measuredFrFlow, 0);
		
	}

}
