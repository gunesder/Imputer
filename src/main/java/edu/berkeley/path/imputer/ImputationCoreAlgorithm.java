package edu.berkeley.path.imputer;

import java.io.IOException;
import java.util.*;

//import jxl.read.biff.BiffException;

import lpsolve.LpSolve;
import lpsolve.LpSolveException;



public class ImputationCoreAlgorithm {
	
	// fields
		// primary fields
		private LinkedList<Cell> cellData = new LinkedList<Cell>();
		private HashMap<Integer,Detector> detectorList = new HashMap<Integer,Detector>();
		private Double demandTimeStep = 5.0/60.0; // [hours] 5 minute default
		private Double simulationTimeStep = 5.0/60.0/60.0; // [hours] 5 second default
		private boolean downBoundaryCongested = false;
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
		private Double[][] measuredDensity = new Double[288][1];
		private Double[][] measuredFlow = new Double[288][1];
		private Double[][] measuredSpeed = new Double[288][1];
		private Double[][] measuredOrFlow = new Double[288][1];
		private Double[][] measuredFrFlow = new Double[288][1];		
		// learning algorithm parameters
	
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
	
	public ImputationCoreAlgorithm(LinkedList<Cell> cells, HashMap<Integer,Detector> detectors, boolean flwBoundary, Double simTimeStep){
		this.cellData = cells;
		this.detectorList = detectors;
		this.downBoundaryCongested = flwBoundary;
		this.simulationTimeStep = simTimeStep;
		this.initializeDataMatrices();
	}
	
	public ImputationCoreAlgorithm(LinkedList<Cell> cells, HashMap<Integer,Detector> detectors, boolean flwBoundary, Double simTimeStep, Double demTimeStep){
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
			Double dummyLength = 0.0;
			for(Link l:c.getLinks()){
				dummyLength = dummyLength + l.getLength();
			}
			linkLength.add(dummyLength);
			
			// Lanes
			lane.add(c.getLinks().getFirst().getLanesML()+c.getLinks().getFirst().getLanesHOV());
			
			// Flow Measurements
			if (separateHOV){
				Double[] tempFLW_ML = c.getDetectorML().getSmoothedDataArray("flow");
				Double[] tempFLW_HOV = c.getDetectorHOV().getSmoothedDataArray("flow");
				Double[] tempFLW = MyUtilities.addVectors(tempFLW_ML,tempFLW_HOV);
				tempFLW = MyUtilities.scaleVector(tempFLW,this.simulationTimeStep*(c.getLinks().getFirst().getLanesML()+c.getLinks().getFirst().getLanesHOV()));
				measuredFlow = MyUtilities.assignColumn(measuredFlow, tempFLW, i);
			} else {
				measuredFlow = MyUtilities.assignColumn(measuredFlow, MyUtilities.scaleVector(c.getDetectorML().getSmoothedDataArray("flow"),this.simulationTimeStep*c.getLinks().getFirst().getLanesML()), i);
			}
			
			// Speed Measurements
			if (separateHOV){
				Double[] tempSPD_ML = c.getDetectorML().getSmoothedDataArray("speed");
				Double[] tempSPD_HOV = c.getDetectorHOV().getSmoothedDataArray("speed");
				tempSPD_ML = MyUtilities.scaleVector(tempSPD_ML,(double) c.getLinks().getFirst().getLanesML());
				tempSPD_HOV = MyUtilities.scaleVector(tempSPD_HOV,(double) c.getLinks().getFirst().getLanesHOV());
				Double[] tempSPD = MyUtilities.addVectors(tempSPD_ML, tempSPD_HOV);
				tempSPD = MyUtilities.scaleVector(tempSPD, (double) (1/(c.getLinks().getFirst().getLanesML()+c.getLinks().getFirst().getLanesHOV())));
				measuredSpeed = MyUtilities.assignColumn(measuredSpeed, tempSPD, i);
			} else {
				measuredSpeed = MyUtilities.assignColumn(measuredSpeed, c.getDetectorML().getSmoothedDataArray("speed"), i);
			}
			
			// Density Measurements
			if (separateHOV){
				Double[] tempDTY_ML = c.getDetectorML().getSmoothedDataArray("density");
				Double[] tempDTY_HOV = c.getDetectorHOV().getSmoothedDataArray("density");
				Double[] tempDTY = MyUtilities.addVectors(tempDTY_ML,tempDTY_HOV);
				measuredDensity = MyUtilities.assignColumn(measuredDensity, MyUtilities.scaleVector(tempDTY,dummyLength*(c.getLinks().getFirst().getLanesML()+c.getLinks().getFirst().getLanesHOV())), i);
			} else {
				measuredDensity = MyUtilities.assignColumn(measuredDensity, MyUtilities.scaleVector(c.getDetectorML().getSmoothedDataArray("density"),dummyLength*c.getLinks().getFirst().getLanesML()), i);
			}
			
			// Fundamental Diagram Parameters
			if (separateHOV) { // case where HOV lane detection is separate
				
				// Mainline
				Double rhocrit_ML = c.getDetectorML().getFdParams().getRho_crit()*dummyLength;
				Double rhojam_ML = c.getDetectorML().getFdParams().getRho_jam()*dummyLength;
				Double qmax_ML = c.getDetectorML().getFdParams().getQ_max()*this.simulationTimeStep;
				
				// HOV
				Double rhocrit_HOV = c.getDetectorHOV().getFdParams().getRho_crit()*dummyLength;
				Double rhojam_HOV = c.getDetectorHOV().getFdParams().getRho_jam()*dummyLength;
				Double qmax_HOV = c.getDetectorHOV().getFdParams().getQ_max()*this.simulationTimeStep;
				
				// Aggregate
				vf.add((qmax_ML+qmax_HOV)/(rhocrit_ML+rhocrit_HOV));
				w.add((qmax_ML+qmax_HOV)/(rhojam_ML+rhojam_HOV-rhocrit_ML-rhocrit_HOV));
				rhojam.add((rhojam_ML+rhojam_HOV)*lane.get(i));
				qmax.add((qmax_ML+qmax_HOV)*lane.get(i));
			
			} else {
				
				vf.add(c.getDetectorML().getFdParams().getVf()*this.simulationTimeStep/dummyLength);
				w.add(c.getDetectorML().getFdParams().getW()*this.simulationTimeStep/dummyLength);
				rhojam.add(c.getDetectorML().getFdParams().getRho_jam()*dummyLength*lane.get(i));
				qmax.add(c.getDetectorML().getFdParams().getQ_max()*this.simulationTimeStep*lane.get(i));
				
			}
			
			// Known Ramp Flows
			Double[] dummyVector = MyUtilities.sumColumns(c.getMeasuredOnrampFlow());
			measuredOrFlow = MyUtilities.assignColumn(measuredOrFlow, dummyVector, i);
			
			Double[] dummyVector2 = MyUtilities.sumColumns(c.getMeasuredOfframpFlow());
			measuredFrFlow = MyUtilities.assignColumn(measuredFrFlow, dummyVector2, i);
			
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
		
		// manually fixing FDparams for testing purposes:
//		qmax.clear();
//		qmax.add(7.6575);
//		qmax.add(7.4796);
//		qmax.add(7.3884);
//		
//		vf.clear();
//		vf.add(0.3158);
//		vf.add(0.1363);
//		vf.add(0.1150);
//		
//		w.clear();
//		w.add(0.0411);
//		w.add(0.0272);
//		w.add(0.0224);
//		
//		rhojam.clear();
//		rhojam.add(210.6119);
//		rhojam.add(330.1817);
//		rhojam.add(393.6067);
		
		// loading and fixing interpolated measurements from text files (exported from matlab)
//		try {
//			measuredSpeed = MyUtilities.read2DArrayFromExcel("C:\\Users\\gsr04\\Workspace\\imputer\\Speed.xls");
//			measuredFlow = MyUtilities.read2DArrayFromExcel("C:\\Users\\gsr04\\Workspace\\imputer\\Flow.xls");
//			measuredDensity = MyUtilities.read2DArrayFromExcel("C:\\Users\\gsr04\\Workspace\\imputer\\Density.xls");
//		} catch (BiffException e1) {
//			// TODO Auto-generated catch block
//			e1.printStackTrace();
//		} catch (IOException e1) {
//			// TODO Auto-generated catch block
//			e1.printStackTrace();
//		}
		
		
		
		// variables needed for the downstream boundary congested case (line 126 of matlab code)
		Double[] boundaryVelocity = MyUtilities.fetchColumn(measuredSpeed, measuredSpeed[0].length-1);
		Double boundaryVF = this.vf.get(this.vf.size()-1)*this.linkLength.get(this.linkLength.size()-1)/this.simulationTimeStep; 
		
		// additional post-processing of variables (line 165 of matlab code)
		Double[] orFLW_save = MyUtilities.fetchColumn(measuredOrFlow, 0);
		Double[] inputFLW = MyUtilities.addVectors(MyUtilities.fetchColumn(measuredFlow, 0),orFLW_save);
		// skipped the assignment in line 169 of matlab
		
		// n'th cell offramp is paired with the n+1'th cell onramp at the node level. the following manipulations are made to take this into account
		imputeOR.remove(0);
		imputeFR.remove(imputeFR.size()-1);
		orPresent.remove(0);
		frPresent.remove(frPresent.size()-1);
		
		Double[][] OrFlow = MyUtilities.removeColumn(measuredOrFlow, 0);
		Double[][] FrFlow = MyUtilities.removeColumn(measuredFrFlow, measuredFrFlow[0].length-1);
		
		Double[][] orFlow_Giv = measuredOrFlow; // probably superfluous but keeping up with the matlab code for now
		Double[][] frFlow_Giv = measuredFrFlow; // probably superfluous but keeping up with the matlab code for now
		
		ArrayList<Boolean> impute = new ArrayList<Boolean>();
		for (int j=0;j<imputeOR.size();j++){
			impute.add(imputeOR.get(j)||imputeFR.get(j));
		}
		
		int numberOfNodes = measuredDensity[0].length + 1;
		
		Double[][] BETA = new Double[FrFlow.length][FrFlow[0].length+1];
		Double[][] dj = new Double[FrFlow.length][FrFlow[0].length+1];
		BETA = MyUtilities.zerosMatrix(FrFlow.length,FrFlow[0].length+1);
		dj = MyUtilities.zerosMatrix(FrFlow.length,FrFlow[0].length+1);
		Double[][] OrINP = new Double[dj.length][dj[0].length];
		Double[] dprev = new Double[numberOfNodes-2];
		Arrays.fill(dprev, 0.0);
		Double[][] ErrA = new Double[BETA.length][BETA[0].length];
		Double[][] ErrB = new Double[BETA.length][BETA[0].length];
		ErrA = MyUtilities.zerosMatrix(BETA.length, BETA[0].length);
		ErrB = MyUtilities.zerosMatrix(BETA.length, BETA[0].length);
		
		// skipped line 198 for now
		
		// *********************************************************************************************
		// ************** Derive dj and BETA from known ramp flows *************************************
		// *********************************************************************************************
		
		Double[] STime = MyUtilities.createIncrementVector(0.0, 24-this.simulationTimeStep, this.simulationTimeStep);
		for (int ii=1;ii<STime.length;ii++){
			
			for (int j=0;j<numberOfNodes-2;j++){
				
				if (!impute.get(j)){
					
					 Double Term1 = this.qmax.get(j+1) < this.w.get(j+1)*(this.rhojam.get(j+1) - measuredDensity[ii][j+1]) ? this.qmax.get(j+1) : this.w.get(j+1)*(this.rhojam.get(j+1)-measuredDensity[ii][j+1]);
					 Double Term2 = this.qmax.get(j) < measuredDensity[ii][j]*this.vf.get(j) ? this.qmax.get(j) : measuredDensity[ii][j]*this.vf.get(j);
					 
					 Double r = OrFlow[ii][j];
					 Double s = FrFlow[ii][j] < Term2 ? FrFlow[ii][j] : Term2;
					 
					 if (Term2+r-s>Term1){
						 
						 Double a1 = Term1 - r;
						 Double a2 = Term2*r;
						 Double c1 = a2;
						 Double b1 = -s;
						 Double b2 = Term2*(Term1+s);
						 Double c2 = s*Term2;
						 int[] colno = new int[4];
				         double[] row = new double[4];
						 					 
						 try {
							 
							LpSolve solver = LpSolve.makeLp(0, 4); // 0 constraints, 4 variables to start with
							solver.strSetObjFn("1 1 0 0"); // cost is the sum of the first two variables
							
							solver.setAddRowmode(true);  /* makes building the model faster if it is done rows by row */

				            /* construct first row (-x1 + a1x3 +a2x4 <= c1) */
				            j = 0;

				            colno[j] = 1; /* first column */
				            row[j++] = -1.0;

				            colno[j] = 2; /* second column */
				            row[j++] = 0.0;
				            
				            colno[j] = 3; /* third column */
				            row[j++] = a1;

				            colno[j] = 4; /* fourth column */
				            row[j++] = a2;

				            /* add the row to lpsolve */
				            solver.addConstraintex(j, row, colno, LpSolve.LE, c1);
				            
				            /* construct second row (-x1 - a1x3 - a2x4 <= -c1) */
				            j = 0;

				            colno[j] = 1; /* first column */
				            row[j++] = -1;

				            colno[j] = 2; /* second column */
				            row[j++] = 0;
				            
				            colno[j] = 3; /* third column */
				            row[j++] = -a1;

				            colno[j] = 4; /* fourth column */
				            row[j++] = -a2;

				            /* add the row to lpsolve */
				            solver.addConstraintex(j, row, colno, LpSolve.LE, -c1);
				            
				            /* construct third row (-x2 + b1x3 + b2x4 <= c2) */
				            j = 0;

				            colno[j] = 1; /* first column */
				            row[j++] = 0;

				            colno[j] = 2; /* second column */
				            row[j++] = -1;
				            
				            colno[j] = 3; /* third column */
				            row[j++] = b1;

				            colno[j] = 4; /* fourth column */
				            row[j++] = b2;

				            /* add the row to lpsolve */
				            solver.addConstraintex(j, row, colno, LpSolve.LE, c2);
				            
				            /* construct fourth row (-x2 - b1x3 - b2x4 <= -c2) */
				            j = 0;

				            colno[j] = 1; /* first column */
				            row[j++] = 0;

				            colno[j] = 2; /* second column */
				            row[j++] = -1;
				            
				            colno[j] = 3; /* third column */
				            row[j++] = -b1;

				            colno[j] = 4; /* fourth column */
				            row[j++] = -b2;

				            /* add the row to lpsolve */
				            solver.addConstraintex(j, row, colno, LpSolve.LE, -c2);
				            
				            /* construct fifth row (x4 <= 1) */
				            j = 0;

				            colno[j] = 1; /* first column */
				            row[j++] = 0;

				            colno[j] = 2; /* second column */
				            row[j++] = 0;
				            
				            colno[j] = 3; /* third column */
				            row[j++] = 0;

				            colno[j] = 4; /* fourth column */
				            row[j++] = 1;

				            /* add the row to lpsolve */
				            solver.addConstraintex(j, row, colno, LpSolve.LE, 1);
				            
				            /* construct sixth row (-x3 <= -max(r,dprev)) */
				            j = 0;

				            colno[j] = 1; /* first column */
				            row[j++] = 0;

				            colno[j] = 2; /* second column */
				            row[j++] = 0;
				            
				            colno[j] = 3; /* third column */
				            row[j++] = -1;

				            colno[j] = 4; /* fourth column */
				            row[j++] = 0;

				            /* add the row to lpsolve */
				            solver.addConstraintex(j, row, colno, LpSolve.LE, (r > dprev[j] ? r : dprev[j]));
				            
				            /* construct seventh row (-x4 <= 0) */
				            j = 0;

				            colno[j] = 1; /* first column */
				            row[j++] = 0;

				            colno[j] = 2; /* second column */
				            row[j++] = 0;
				            
				            colno[j] = 3; /* third column */
				            row[j++] = 0;

				            colno[j] = 4; /* fourth column */
				            row[j++] = -1;

				            /* add the row to lpsolve */
				            solver.addConstraintex(j, row, colno, LpSolve.LE, 0);
				            
				            solver.setAddRowmode(false); /* rowmode should be turned off again when done building the model */
				            
				            /* Now let lpsolve calculate a solution */
				            int ret = solver.solve();
				            double[] optVariables = new double[4];
				            solver.getVariables(optVariables);
				            
				            /* Assign optimal variables to corresponding local variables */
				            dj[ii][j] = optVariables[2];
				            BETA[ii][j] = optVariables[3];
				            ErrA[ii][j] = optVariables[0];
				            ErrB[ii][j] = optVariables[1];
				            OrINP[ii][j] = dj[ii][j] - dprev[j];
				            dprev[j] = dj[ii][j] - dj[ii][j]*Term1/(Term2*(1-optVariables[3])+dj[ii][j]);
							
						} catch (LpSolveException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						 
					 } else {
						dj[ii][j] = r;
						dprev[j] = 0.0;
						OrINP[ii][j] = dj[ii][j] - dprev[j];
						BETA[ii][j] = s/(0.0000000001 > Term2 ? 0.0000000001 : Term2);
					 }
					 
				}
			}
		}
		
		Double[][] Demand = dj.clone(); 
		
		// *********************************************************************************************
		// ************** Check that ramps are present whenever imputation is enabled ******************
		// *********************************************************************************************
		
		ArrayList<Double> lowerBounds = new ArrayList<Double>(); for (int j=0;j<numberOfNodes-2;j++) lowerBounds.add(0.0);
		ArrayList<Double> upperBounds = new ArrayList<Double>(); for (int j=0;j<numberOfNodes-2;j++) upperBounds.add(0.0);
		for (int ind = 0;ind < numberOfNodes - 2; ind++){
			if (frPresent.get(ind)) lowerBounds.set(ind, 0.0); else lowerBounds.set(ind, 1.0);
			if (orPresent.get(ind)) upperBounds.set(ind, 0.0); else upperBounds.set(ind, 1.0);
			if (!frPresent.get(ind) && !orPresent.get(ind) && impute.get(ind)){
				// some error message given here originally
			}			
		}
		if (downBoundaryCongested){
			lowerBounds.add(1.0);
			upperBounds.add(0.0);
		}
		
		// *********************************************************************************************
		// ************************* Learning Algorithm ************************************************
		// *********************************************************************************************
		
		/* Learning algorithm front matter, gains, user-defined settings, etc. (hard coded for now) */
		Double GM = 40.0; // these are user-defined gains for the adaptive learning algorithm and are not necessarily all very sensitive
		Double G1 = 1*GM;
		Double G2 = 0.001*GM;
		Double percTol = 0.001;
		Double percTol2 = 0.2;
		int maxLim = 100;
		int iterMax = 25;
		int[] iterTrigger = {4,8,11,15,19}; // iteration indices when the trigger algorithm kicks in
		Double derivativeBound = 2.0;
		Double startBound = 0.0;
		
		// initialize c: this is the effective demand vector into each cell
		Double[][] c = new Double[STime.length][cellData.size()-1];
		if (downBoundaryCongested){
			c = new Double[STime.length][cellData.size()];
		}
		Double[] dummy = new Double[measuredDensity.length];
		for (int j=1;j<numberOfNodes-1;j++){
			dummy = MyUtilities.scaleVector(MyUtilities.onesVector(measuredDensity.length), rhojam.get(j));
			dummy = MyUtilities.addVectors(dummy,MyUtilities.scaleVector(MyUtilities.fetchColumn(measuredDensity, j),-1.0));
			c = MyUtilities.assignColumn(c, MyUtilities.scaleVector(dummy, 1.5*w.get(j)), j-1);
		}
		
		if (downBoundaryCongested){
			MyUtilities.assignColumn(c, MyUtilities.scaleVector(MyUtilities.onesVector(c.length), qmax.get(numberOfNodes-2)), numberOfNodes-2);
			impute.add(true);
		}
		
		Double[][] cBest = c;
		int iterBound_D_Beta = iterMax+1;
		int startIterBound = 10;
		
		boolean[] flags = new boolean[numberOfNodes];
		
		// initialize csave: this is the effective demand vector arraylist for keeping track of the c matrix between iterations
		LinkedList<Double[][]> csave = new LinkedList<Double[][]>();
		LinkedList<Double> MError = new LinkedList<Double>(); 
		Double[][] Nh = MyUtilities.zerosMatrix(STime.length, cellData.size());
		
		// Learning Algorithm Loop
		for (int iter = 0;iter < iterMax;iter++){
			
			flags[0] = true;
			
			Double[][] mode = MyUtilities.zerosMatrix(c.length, numberOfNodes-1);
			
			Double InQ = 0.0;
			Nh = MyUtilities.assignRow(Nh, measuredDensity[0], 0);
			Double[] cj = c[0];
			Double[] Limit = new Double[cellData.size()];
			
			for (int k=0;k<STime.length-1;k++){ // starts at line 328 and ends at line 584
				
				for (int j=0;j<Limit.length;j++){
					Limit[j] = qmax.get(j) < w.get(j)*(rhojam.get(j)-Nh[k][j]) ? qmax.get(j) : w.get(j)*(rhojam.get(j)-Nh[k][j]);
				}
				
				Double[] cjprev = cj.clone();
				cj = c[k];
				
				/* This whole section does what line 333 does in Matlab, may need to make into a utility method if it is needed somewhere else*/
				Double[] dummy1 = MyUtilities.subtractVectors(cjprev,MyUtilities.scaleVector(MyUtilities.onesVector(cjprev.length), derivativeBound));
				Double[] dummy2 = MyUtilities.addVectors(cjprev, MyUtilities.scaleVector(MyUtilities.onesVector(cjprev.length), derivativeBound));
				Double[] dummy3 = new Double[cj.length];
				
				for (i=0;i<cj.length;i++){
					dummy3[i] = cj[i] < dummy2[i] ? cj[i] : dummy2[i];
				}
				
				for (i=0;i<cj.length;i++){
					cj[i] = dummy1[i] > dummy3[i] ? dummy1[i] : dummy3[i];
				}
				/* This whole section does what line 333 does in Matlab, may need to make into a utility method if it is needed somewhere else*/
				
				for (int j=0;j<numberOfNodes-2;j++){
					
					Double NjVfj = qmax.get(j) < Nh[k][j]*vf.get(j) ? qmax.get(j) : Nh[k][j]*vf.get(j);
					
					if(!impute.get(j)){
						cj[j] = NjVfj*(1-BETA[k][j]) + Demand[k][j];
					} else {
						if (iter>startIterBound){
							if (lowerBounds.get(j) !=0){
								cj[j] = cj[j] > NjVfj ? cj[j] : NjVfj;
							}
							if (upperBounds.get(j) != 0){
								cj[j] = cj[j] < NjVfj ? cj[j] : NjVfj;
							}
						}
					}
				}
				
				for (int j=1;j<numberOfNodes-1;j++){
					if (impute.get(j-1)){
						flags[j] = cj[j-1]>=Limit[j]*(1-percTol);
					} else {
						flags[j] = cj[j-1]>=Limit[j];
					}
				}
				flags[numberOfNodes-1] = false;
				
				InQ = InQ + inputFLW[k];
				Double inFlow = Limit[0] < InQ ? Limit[0] : InQ;
				InQ = InQ - inFlow;
				
				boolean needBoundaryImpute = downBoundaryCongested && (boundaryVelocity[k]<boundaryVF*0.9);
				if (needBoundaryImpute){
					flags[numberOfNodes-1] = true;
				}
				
				int j=-1;
				while (j<numberOfNodes-2){ // starts at line 374 and ends at line 580
					j++;
					Double NjVfj = Nh[k][j]*vf.get(j) < qmax.get(j) ? Nh[k][j]*vf.get(j) : qmax.get(j);
					Double Njm1Vfjm1 = 0.0; // never used, just a workaround. It is assigned a new value for each node except for the first one, which doesn't undergo imputation anyway
					if (j>=1){
						Njm1Vfjm1 = qmax.get(j-1) < Nh[k][j-1]*vf.get(j-1) ? qmax.get(j-1) : Nh[k][j-1]*vf.get(j-1);	
					}
					Double Wj = qmax.get(j) < w.get(j)*(rhojam.get(j)-Nh[k][j]) ? qmax.get(j) : w.get(j)*(rhojam.get(j)-Nh[k][j]);
					Double Wjp1 = 0.0;
					if (j == numberOfNodes-2){
						Wjp1 = qmax.get(j);
					} else {
						Wjp1 = qmax.get(j+1) < w.get(j+1)*(rhojam.get(j+1)-Nh[k][j+1]) ? qmax.get(j+1) : w.get(j+1)*(rhojam.get(j+1)-Nh[k][j+1]);
					}
					
					if (!flags[j] && !flags[j+1]){ // starts at line 391 and ends at line 577
						mode[k][j] = 1.0;
						
						Double BoundL = lowerBounds.get(j-1);
						Double BoundU = upperBounds.get(j-1);
						
						Double Nh0 = measuredDensity[k+1][j] - (Nh[k][j]-NjVfj+cj[j-1]);
						Double NHt = Nh0/(1+G1);
						
						if (impute.get(j-1)) { // Line 402 to 424
							
							cj[j-1] = cj[j-1] + G1*NHt;
							
							Double UBound = 0.0;
							Double LBound = 0.0;
							if (iter > iterBound_D_Beta){
								UBound = Wj < Njm1Vfjm1 + maxD.get(j-1) ? Wj : Njm1Vfjm1 + maxD.get(j-1);
								LBound = 0.01 > (Wj < Njm1Vfjm1*(1-maxBeta.get(j-1)) ? Wj : Njm1Vfjm1*(1-maxBeta.get(j-1))) ? 0.01 : (Wj < Njm1Vfjm1*(1-maxBeta.get(j-1)) ? Wj : Njm1Vfjm1*(1-maxBeta.get(j-1)));
							} else {
								UBound = Wj;
								LBound = 0.01;
							}
							
							if (BoundU == 1){
								UBound = UBound < Njm1Vfjm1 ? UBound : Njm1Vfjm1;
							} else if (BoundL == 1){
								LBound = LBound > Njm1Vfjm1 ? LBound : Njm1Vfjm1;
							}
							
							if (k>startBound){
								LBound = LBound > cjprev[j-1]-derivativeBound ? LBound : cjprev[j-1]-derivativeBound;
								UBound = UBound < cjprev[j-1]+derivativeBound ? UBound : cjprev[j-1]+derivativeBound;
							}
							
							cj[j-1] = (cj[j-1] < UBound ? cj[j-1] : UBound) > LBound ? (cj[j-1] < UBound ? cj[j-1] : UBound) : LBound;
							
						}
						
						Nh[k+1][j] = Nh[k][j] - NjVfj + cj[j-1];
						NHt = measuredDensity[k+1][j] - Nh[k+1][j];
						
					} else if (!flags[j] && flags[j+1]){
						mode[k][j] = 2.0;
						
						if (impute.get(j) && cj[j] > Wjp1*(1+percTol) && cj[j] < Wjp1*(1+percTol) && (Nh[k][j]+cj[j-1]-NjVfj)>measuredDensity[k+1][j]){
							cj[j] = cj[j] < Wjp1*(1-percTol/100) ? cj[j] : Wjp1*(1-percTol/100);
							flags[j+1] = false;
							j = j-1;
							continue;							
						}
						
						cj[j] = cj[j] > Wjp1 ? cj[j] : Wjp1;
						
						Double BoundL1 = lowerBounds.get(j-1);
						Double BoundU1 = upperBounds.get(j-1);
						Double BoundL2 = lowerBounds.get(j);
						Double BoundU2 = upperBounds.get(j);
						
						Double Nh0 = measuredDensity[k+1][j] - (Nh[k][j]+cj[j-1]-Wjp1*NjVfj/cj[j]);
						Double NHt = Nh0 / (1+G1*(impute.get(j-1) ? 1 : 0)+G2*Wjp1*NjVfj*(impute.get(j) ? 1 : 0));
						
						if (impute.get(j)){
							cj[j] = 1 / (1/cj[j]-(G2*NHt < 1/cj[j]*0.99 ? G2*NHt : 1/cj[j]*0.99));
							
							Double UBound2 = 0.0;
							if (BoundU2 != 0){
								UBound2 =  (maxLim*Wjp1 < NjVfj ? maxLim*Wjp1 : NjVfj) > 0.01 ? (maxLim*Wjp1 < NjVfj ? maxLim*Wjp1 : NjVfj) : 0.01;
							} else if (iter > iterBound_D_Beta){
								UBound2 = maxLim*Wjp1 < (Wjp1 > NjVfj+maxD.get(j) ? Wjp1 : NjVfj+maxD.get(j)) ? maxLim*Wjp1 : (Wjp1 > NjVfj+maxD.get(j) ? Wjp1 : NjVfj+maxD.get(j));
							} else {
								UBound2 = maxLim*Wjp1;
							}
							
							Double LBound2 = NjVfj*BoundL2 > Wjp1 ? NjVfj*BoundL2 : Wjp1;
							
							if (iter > iterBound_D_Beta){
								LBound2 = LBound2 > NjVfj*(1-maxBeta.get(j)) ? LBound2 : NjVfj*(1-maxBeta.get(j));
							} else {
								LBound2 = LBound2 > 0.01 ? LBound2 : 0.01; 
							}
							
							if (k > startBound){
								LBound2 = LBound2 > cjprev[j]-derivativeBound ? LBound2 : cjprev[j]-derivativeBound;
								UBound2 = UBound2 < cjprev[j]+derivativeBound ? UBound2 : cjprev[j]+derivativeBound;
							}
							
							cj[j] = (cj[j] < UBound2 ? cj[j] : UBound2) > LBound2 ? (cj[j] < UBound2 ? cj[j] : UBound2) : LBound2;							
							
						}
						
						if (impute.get(j-1)){
							cj[j-1] = cj[j-1] + G1*NHt;
							
							Double UBound1 = 0.0;
							Double LBound1 = 0.0;
							if (iter > iterBound_D_Beta){
								UBound1 = Wj < Njm1Vfjm1+maxD.get(j-1) ? Wj : Njm1Vfjm1+maxD.get(j-1);
								LBound1 = 0.01 > (Wj < Njm1Vfjm1*(1-maxBeta.get(j-1)) ? Wj : Njm1Vfjm1*(1-maxBeta.get(j-1))) ? 0.01 : (Wj < Njm1Vfjm1*(1-maxBeta.get(j-1)) ? Wj : Njm1Vfjm1*(1-maxBeta.get(j-1)));
							} else {
								UBound1 = Wj;
								LBound1 = 0.01;
							}
							
							if (BoundU1 == 1){
								UBound1 = Njm1Vfjm1;
							} else if (BoundL1 == 1){
								LBound1 = LBound1 > Njm1Vfjm1 ? LBound1 : Njm1Vfjm1;
							}
							
							if (k > startBound){
								LBound1 = LBound1 > cjprev[j-1] - derivativeBound ? LBound1 : cjprev[j-1] - derivativeBound;
								UBound1 = UBound1 < cjprev[j-1] + derivativeBound ? UBound1 : cjprev[j-1] + derivativeBound;
							}
							
							cj[j] = (cj[j-1] < UBound1 ? cj[j-1] : UBound1) > LBound1 ? (cj[j-1] < UBound1 ? cj[j-1] : UBound1) : LBound1;
								
						}
						
						Nh[k+1][j] = Nh[k][j] + cj[j-1] - Wjp1*NjVfj/cj[j];
						NHt = measuredDensity[k+1][j]-Nh[k+1][j];						
						
					} else if (flags[j] && !flags[j+1]){
						
						mode[k][j] = 3.0;
						if (j==0){
							Nh[k+1][j] = Nh[k][j] + inFlow - (Nh[k][j]*vf.get(j) < qmax.get(j) ? Nh[k][j]*vf.get(j) : qmax.get(j));
						} else {
							Nh[k+1][j] = Nh[k][j] + (qmax.get(j) < w.get(j)*(rhojam.get(j)-Nh[k][j]) ? qmax.get(j) : w.get(j)*(rhojam.get(j)-Nh[k][j])) - (Nh[k][j]*vf.get(j) < qmax.get(j) ? Nh[k][j]*vf.get(j) : qmax.get(j));
						}
						
					} else if (flags[j] && flags[j+1]){
						mode[k][j] = 4.0;
						
						if (impute.get(j)){
							Double BoundL = lowerBounds.get(j);
							Double BoundU = upperBounds.get(j);
							
							Double Nh0 = 0.0;
							if (j==0){
								Nh0 = measuredDensity[k+1][j] - (Nh[k][j]+inFlow-Wjp1*NjVfj/cj[j]);
							} else {
								Nh0 = measuredDensity[k+1][j] - (Nh[k][j]+Wj-Wjp1*NjVfj/cj[j]);
							}
							
							if (cj[j]>Wjp1*(1-percTol) && cj[j]<Wjp1*(1+percTol) && Nh0<0){
								cj[j] = cj[j] < Wjp1*(1-percTol/100) ? cj[j] : Wjp1*(1-percTol/100);
								flags[j+1] = false;
								j = j-1;
								continue;							
							}
							
							cj[j] = cj[j] > Wjp1 ? cj[j] : Wjp1;
							
							Double NHt = Nh0/(1+G2*Wjp1*NjVfj);
							
							cj[j] = 1/(1/cj[j]-(G2*NHt < 1/cj[j]*0.99 ? G2*NHt : 1/cj[j]*0.99));
							
							Double UBound = 0.0;
							if (BoundU != 0){
								UBound = (maxLim*Wjp1 < NjVfj ? maxLim*Wjp1 : NjVfj) > 0.01 ? (maxLim*Wjp1 < NjVfj ? maxLim*Wjp1 : NjVfj) : 0.01;
							} else if (iter > iterBound_D_Beta){
								UBound = maxLim*Wjp1 < (Wjp1 > NjVfj*maxD.get(j) ? Wjp1 : NjVfj*maxD.get(j)) ? maxLim*Wjp1 : (Wjp1 > NjVfj*maxD.get(j) ? Wjp1 : NjVfj*maxD.get(j));
							} else {
								UBound = maxLim*Wjp1;
							}
							
							Double LBound = NjVfj*BoundL > Wjp1 ? NjVfj*BoundL : Wjp1;
							if (iter > iterBound_D_Beta){
								LBound = LBound > NjVfj*(1-maxBeta.get(j)) ? LBound : NjVfj*(1-maxBeta.get(j));
							} else {
								LBound = LBound > 0.01 ? LBound : 0.01;
							}
							
							if (k > startBound){
								LBound = LBound > cjprev[j] - derivativeBound ? LBound : cjprev[j] - derivativeBound;
								UBound = UBound < cjprev[j] + derivativeBound ? UBound : cjprev[j] + derivativeBound;
							}
							
							cj[j] = (cj[j] < UBound ? cj[j] : UBound) > LBound ? (cj[j] < UBound ? cj[j] : UBound) : LBound;
							
						}
						
						if (j==0){
							Nh[k+1][j] = Nh[k][j] + inFlow - Wjp1*NjVfj/cj[j];
						} else {
							Nh[k+1][j] = Nh[k][j] + Wj - Wjp1*NjVfj/cj[j];
						}
						
						Double NHt = measuredDensity[k+1][j] - Nh[k+1][j];
						
					}
					
				} // end of while loop over nodes
				
				for (i=0;i<cj.length;i++){
					cj[i] = 0.01 > cj[i] ? 0.01 : cj[i];
				}
				
				MyUtilities.assignRow(c, cj, k);
				
			} // end of for loop over time
			
			csave.add(iter, c);
			// skip the check in lines 587 to 589
			Double[][] dummy1 = MyUtilities.addMatrices(measuredDensity, MyUtilities.scaleMatrix(Nh, -1.0));
			dummy1 = MyUtilities.matrixAbsValue(dummy1);
			Double[] dummy2 = MyUtilities.meanColumns(dummy1);
			Double[] dummy3 = MyUtilities.meanColumns(measuredDensity);
			MError.add(iter, MyUtilities.meanVector(dummy2)/MyUtilities.meanVector(dummy3));
			
			System.out.println("Density Error: " + MError.get(iter)*100 + ", iteration: " + (iter+1) + " of " + iterMax);
			
			boolean flag = true;
			for (int j=0;j<=iter;j++){
				if (MError.get(iter)>MError.get(j)){
					flag = false;
				}
			}
			if (iter == startIterBound || (iter>startIterBound && flag )){
				cBest = c;
			}
			
			// Trigger algorithm to be triggered at pre-determined iterations
			boolean flag2 = false;
			for (int j=0;j<iterTrigger.length;j++){
				if (iter == iterTrigger[j]){
					flag2 = true;
				}
			}
			
			if (flag2){ // Trigger Algorithm block
				
				// Find regions that have bottlenecks (i.e. mode 3) Term I line 602
				boolean[][] ModeFilter = new boolean[mode.length][mode[0].length];
				for (int row=0;row<ModeFilter.length;row++){ 
					for (int col=0;col<ModeFilter[0].length;col++){
						if (mode[row][col] == 3){
							ModeFilter[row][col] = true; 
						} else {
							ModeFilter[row][col] = false;
						}
					}
				}
				
				if (downBoundaryCongested){
					// Breaking Down line 604 into parts (Term I is the same as ModeFilter above)
					// Term II, Part 1 (Mode == 1)
					boolean[][] term2part1 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term2part1.length;row++){
						for (int col=0;col<term2part1[0].length;col++){
							if (mode[row][col] == 1){
								term2part1[row][col] = true; 
							} else {
								term2part1[row][col] = false;
							}
						}
					}
					
					// Term II, Part 2 (ones(size(Mode,1),1)*[0 Impute(1:end-1)]==0)
					ArrayList<Boolean> dummyImpute = new ArrayList<Boolean>();
					dummyImpute = (ArrayList<Boolean>) impute.clone();
					dummyImpute.remove(dummyImpute.size()-1);
					dummyImpute.add(0, false);
					
					Boolean[] dummyImputeArray = dummyImpute.toArray(new Boolean[dummyImpute.size()]);
					boolean[][] term2part2 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term2part2.length;row++){ 
						for (int col=0;col<term2part2[0].length;col++){
							term2part2[row][col] = true && dummyImputeArray[col]; 
						}
					}
					
					for (int row=0;row<term2part2.length;row++){
						for (int col=0;col<term2part1[0].length;col++){
							if (term2part2[row][col] == false){
								term2part2[row][col] = true; 
							} else {
								term2part2[row][col] = false;
							}
						}
					}
					
					// Term II, Part 3 (Density-Nh>0)
					boolean[][] term2part3 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term2part3.length;row++){
						for (int col=0;col<term2part3[0].length;col++){
							if (measuredDensity[row][col] > Nh[row][col]){
								term2part3[row][col] = true; 
							} else {
								term2part3[row][col] = false;
							}
						}
					}
					
					// Term II, all parts together
					boolean[][] term2 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term2.length;row++){
						for (int col=0;col<term2[0].length;col++){
							term2[row][col] = term2part1[row][col] && term2part2[row][col] && term2part3[row][col];
						}
					}
					
					// Term III, Part 1
					boolean[][] term3part1 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term3part1.length;row++){
						for (int col=0;col<term3part1[0].length;col++){
							if (mode[row][col] == 4){
								term3part1[row][col] = true; 
							} else {
								term3part1[row][col] = false;
							}
						}
					}
					
					// Term III, Part 2
					ArrayList<Boolean> dummyImpute2 = new ArrayList<Boolean>();
					dummyImpute2 = (ArrayList<Boolean>) impute.clone();
					dummyImpute2.remove(dummyImpute2.size()-1);
					dummyImpute2.add(dummyImpute2.size()-1, false);
					Boolean[] dummyImputeArray2 = dummyImpute2.toArray(new Boolean[dummyImpute2.size()]);
					boolean[][] term3part2 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term3part2.length;row++){ 
						for (int col=0;col<term3part2[0].length;col++){
							term3part2[row][col] = true && dummyImputeArray2[col]; 
						}
					}
					
					for (int row=0;row<term3part2.length;row++){
						for (int col=0;col<term2part1[0].length;col++){
							if (term3part2[row][col] == false){
								term3part2[row][col] = true; 
							} else {
								term3part2[row][col] = false;
							}
						}
					}
					
					// Term III, Part 3 (Density-Nh<0)
					boolean[][] term3part3 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term3part3.length;row++){
						for (int col=0;col<term3part3[0].length;col++){
							if (measuredDensity[row][col] < Nh[row][col]){
								term3part3[row][col] = true; 
							} else {
								term3part3[row][col] = false;
							}
						}
					}
					
					// Term III, all parts together
					boolean[][] term3 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term3.length;row++){
						for (int col=0;col<term3[0].length;col++){
							term3[row][col] = term3part1[row][col] && term3part2[row][col] && term3part3[row][col];
						}
					}
					
					// All Terms Combined (line 604 Finally)
					for (int row=0;row<ModeFilter.length;row++){ 
						for (int col=0;col<ModeFilter[0].length;col++){
							ModeFilter[row][col] = ModeFilter[row][col] || term2[row][col] || term3[row][col];
						}
					}
				
				} else {
					
					// Term I again is ModeFilter above.
					
					// Term II, Part 1 (Mode == 1)
					boolean[][] term2part1 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term2part1.length;row++){
						for (int col=0;col<term2part1[0].length;col++){
							if (mode[row][col] == 1){
								term2part1[row][col] = true; 
							} else {
								term2part1[row][col] = false;
							}
						}
					}
					
					// Term II, Part 2 (ones(size(Mode,1),1)*[0 Impute(1:end-1)]==0)
					ArrayList<Boolean> dummyImpute = new ArrayList<Boolean>();
					dummyImpute = (ArrayList<Boolean>) impute.clone();
					dummyImpute.add(0, false);
					Boolean[] dummyImputeArray = dummyImpute.toArray(new Boolean[dummyImpute.size()]);
					boolean[][] term2part2 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term2part2.length;row++){ 
						for (int col=0;col<term2part2[0].length;col++){
							term2part2[row][col] = true && dummyImputeArray[col]; 
						}
					}
					
					for (int row=0;row<term2part2.length;row++){
						for (int col=0;col<term2part1[0].length;col++){
							if (term2part2[row][col] == false){
								term2part2[row][col] = true; 
							} else {
								term2part2[row][col] = false;
							}
						}
					}
					
					// Term II, Part 3 (Density-Nh>0)
					boolean[][] term2part3 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term2part3.length;row++){
						for (int col=0;col<term2part3[0].length;col++){
							if (measuredDensity[row][col] > Nh[row][col]){
								term2part3[row][col] = true; 
							} else {
								term2part3[row][col] = false;
							}
						}
					}
					
					// Term II, all parts together
					boolean[][] term2 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term2.length;row++){
						for (int col=0;col<term2[0].length;col++){
							term2[row][col] = term2part1[row][col] && term2part2[row][col] && term2part3[row][col];
						}
					}
					
					// Term III, Part 1
					boolean[][] term3part1 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term3part1.length;row++){
						for (int col=0;col<term3part1[0].length;col++){
							if (mode[row][col] == 4){
								term3part1[row][col] = true; 
							} else {
								term3part1[row][col] = false;
							}
						}
					}
					
					// Term III, Part 2
					ArrayList<Boolean> dummyImpute2 = new ArrayList<Boolean>();
					dummyImpute2 = (ArrayList<Boolean>) impute.clone();
					dummyImpute2.add(dummyImpute2.size()-1, false);
					Boolean[] dummyImputeArray2 = dummyImpute2.toArray(new Boolean[dummyImpute2.size()]);
					boolean[][] term3part2 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term3part2.length;row++){ 
						for (int col=0;col<term3part2[0].length;col++){
							term3part2[row][col] = true && dummyImputeArray2[col]; 
						}
					}
					
					for (int row=0;row<term3part2.length;row++){
						for (int col=0;col<term2part1[0].length;col++){
							if (term3part2[row][col] == false){
								term3part2[row][col] = true; 
							} else {
								term3part2[row][col] = false;
							}
						}
					}
					
					// Term III, Part 3 (Density-Nh<0)
					boolean[][] term3part3 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term3part3.length;row++){
						for (int col=0;col<term3part3[0].length;col++){
							if (measuredDensity[row][col] < Nh[row][col]){
								term3part3[row][col] = true; 
							} else {
								term3part3[row][col] = false;
							}
						}
					}
					
					// Term III, all parts together
					boolean[][] term3 = new boolean[mode.length][mode[0].length];
					for (int row=0;row<term3.length;row++){
						for (int col=0;col<term3[0].length;col++){
							term3[row][col] = term3part1[row][col] && term3part2[row][col] && term3part3[row][col];
						}
					}
					
					// All Terms Combined (line 606 Finally)
					for (int row=0;row<ModeFilter.length;row++){ 
						for (int col=0;col<ModeFilter[0].length;col++){
							ModeFilter[row][col] = ModeFilter[row][col] || term2[row][col] || term3[row][col];
						}
					}
										
				} // line 607
				
				Double[][] BottleError = new Double[ModeFilter.length][ModeFilter[0].length];
				for (int row=0;row<BottleError.length;row++){
					for (int col=0;col<BottleError[0].length;col++){
						if (ModeFilter[row][col] && measuredDensity[row][col]>0.001){
							BottleError[row][col] = (measuredDensity[row][col]-Nh[row][col])/(measuredDensity[row][col] > 0.001 ? measuredDensity[row][col] : 0.001);
						} else {
							BottleError[row][col] = 0.0;
						}
					}
				}
				
				boolean[][] Filter = new boolean[ModeFilter.length][ModeFilter[0].length];
				for (int row=0;row<Filter.length;row++){ 
					for (int col=0;col<Filter[0].length;col++){
						if (col == Filter[0].length - 1){
							Filter[row][col] = false;
						} else {
							Filter[row][col] = true;
						}
					}
				}
				
				// Find the cells and instances where BottleError is greater than 0.05
				ArrayList<Integer> rowIndeces = new ArrayList<Integer>();
				ArrayList<Integer> colIndeces = new ArrayList<Integer>();
				for (int row=0;row<ModeFilter.length;row++){ 
					for (int col=0;col<ModeFilter[0].length;col++){
						if (BottleError[row][col]>0.05 && Filter[row][col]){
							rowIndeces.add(row);
							colIndeces.add(col);
						}
					}
				}
				
				Double limitNew = 0.0;
				for (int row=0;row<rowIndeces.size();row++){
					if (impute.get(colIndeces.get(row))){
						limitNew = Math.min(qmax.get(colIndeces.get(row)+1), w.get(colIndeces.get(row)+1)*(rhojam.get(colIndeces.get(row)+1) - Nh[rowIndeces.get(row)][colIndeces.get(row)+1]))*(1+percTol2);
						c[rowIndeces.get(row)][colIndeces.get(row)] = limitNew;
					}
				}
				
				ArrayList<Integer> rowIndeces2 = new ArrayList<Integer>();
				ArrayList<Integer> colIndeces2 = new ArrayList<Integer>();
				for (int row=0;row<ModeFilter.length;row++){ 
					for (int col=0;col<ModeFilter[0].length;col++){
						if (BottleError[row][col]<-0.05 && Filter[row][col]){
							rowIndeces2.add(row);
							colIndeces2.add(col);
						}
					}
				}
				
				for (int row=0;row<rowIndeces2.size();row++){
					if (colIndeces2.get(row)>0 && impute.get(colIndeces2.get(row)-1)){
						limitNew = Math.min(qmax.get(colIndeces2.get(row)), w.get(colIndeces2.get(row))*(rhojam.get(colIndeces2.get(row)) - Nh[rowIndeces2.get(row)][colIndeces2.get(row)]))*(1-percTol2);
						c[rowIndeces2.get(row)][colIndeces2.get(row)-1] = limitNew;
					}
				}
				
			} // line 634
			
		} // end of for loop over iterations (line 636)
		
		c = cBest;
		
		// *********************************************************************************************
		// ******* Determine the Parameters dj and Beta from the effective demands *********************
		// *********************************************************************************************
		
		boolean Override = false;
		Double[][] Nha = Nh;
		Double[][] BETA1 = MyUtilities.zerosMatrix(STime.length, cellData.size());
		Double[][] dj1 = MyUtilities.zerosMatrix(STime.length, cellData.size());
		Double[][] OrInp = MyUtilities.zerosMatrix(STime.length, cellData.size());
		Double[][] OnrampFlw = MyUtilities.zerosMatrix(STime.length, cellData.size());
		Double[][] FlowBet = measuredFlow;
		FlowBet = MyUtilities.removeColumn(FlowBet, 1);
		
		for (int ii=1;ii<STime.length;ii++){   // line 655 to 758
			Arrays.fill(dprev, 0.0);
			for (int j=0;j<numberOfNodes-2;j++){
				
				if (impute.get(j)){
					
					 Double Term1 = this.qmax.get(j+1) < this.w.get(j+1)*(this.rhojam.get(j+1) - Nha[ii][j+1]) ? this.qmax.get(j+1) : this.w.get(j+1)*(this.rhojam.get(j+1)-measuredDensity[ii][j+1]);
					 Double Term2 = this.qmax.get(j) < Nha[ii][j]*this.vf.get(j) ? this.qmax.get(j) : Nha[ii][j]*this.vf.get(j);
					 
					 if (!orPresent.get(j)){
						 BETA1[ii][j] = Math.max(0, 1-c[ii][j]/Term2);
						 OnrampFlw[ii][j] = 0.0;
						 continue;
					 }
					 
					 Double MultFac = 0.0;
					 
					 if (c[ii][j]>Term1){
						 
						 if (!frPresent.get(j)){
							 dj1[ii][j] = c[ii][j] - Term2;
							 OrInp[ii][j] = dj1[ii][j] - dprev[j];
							 OnrampFlw[ii][j] = dj1[ii][j]*Term1/c[ii][j];
							 dprev[j] = dj1[ii][j] - OnrampFlw[ii][j];
							 continue;
						 }
						 
						 Term1 = Math.min(this.qmax.get(j+1), this.w.get(j+1)*(this.rhojam.get(j+1) - Nha[ii][j+1]));
						 Term2 = Math.min(this.qmax.get(j), Nha[ii][j]*this.vf.get(j));
						 Double fac = Term1/c[ii][j];
						 
						 // Explicit Solution
						 Double d_soln = 0.0; Double fr_soln = 0.0;
						 if (imputeOR.get(j) && imputeFR.get(j) || Override){
							 
							 Double dummyMax = Math.max(0, dprev[j]);
							 dummyMax = Math.max(dummyMax, Term1 - FlowBet[ii][j]);
							 dummyMax = Math.max(dummyMax, c[ii][j]-Term2);
							 d_soln = Math.min(c[ii][j], dummyMax);
							 fr_soln = Term2 + d_soln - c[ii][j];
							 
						 } else if (!imputeOR.get(j)){
							 
							 MultFac = Math.abs(FlowBet[ii][j]/Math.max(100*this.simulationTimeStep, orFlow_Giv[ii][j]))/2;
							 Double d_vertex = 0.0;
							 if (MultFac<1){
								 d_vertex = (Term1 - FlowBet[ii][j])/fac;
							 } else {
								 d_vertex = orFlow_Giv[ii][j]/fac;
							 }
							 
							 Double dummyMax = Math.max(0, dprev[j]);
							 dummyMax = Math.max(dummyMax, d_vertex);
							 dummyMax = Math.max(dummyMax, c[ii][j]-Term2);
							 d_soln = Math.min(c[ii][j], dummyMax);
							 fr_soln = Term2 + d_soln - c[ii][j];							 
							 
						 } else {
							 
							 MultFac = Math.abs(FlowBet[ii][j]/Math.max(100*this.simulationTimeStep, orFlow_Giv[ii][j]))/2;
							 Double d_vertex = 0.0;
							 if (MultFac<1){
								 d_vertex = (Term1 - FlowBet[ii][j])/fac;
							 } else {
								 d_vertex = (Term1 + frFlow_Giv[ii][j] - fac*Term2)/fac;
							 }
							 
							 Double dummyMax = Math.max(0, dprev[j]);
							 dummyMax = Math.max(dummyMax, d_vertex);
							 dummyMax = Math.max(dummyMax, c[ii][j]-Term2);
							 d_soln = Math.min(c[ii][j], dummyMax);
							 fr_soln = Term2 + d_soln - c[ii][j];
							 
						 }
						 
						 dj1[ii][j] = d_soln;
						 BETA1[ii][j] = fr_soln/Term2;
						 OrInp[ii][j] = d_soln-dprev[j];
						 dprev[j] = d_soln - Term1/c[ii][j]*dj1[ii][j];
						 
					 } else {
						 
						 if (!frPresent.get(j)){
							 dj1[ii][j] = c[ii][j] - Term2;
							 OrInp[ii][j] = dj1[ii][j] - dprev[j];
							 OnrampFlw[ii][j] = dj1[ii][j];
							 dprev[j] = 0.0;
							 continue;
						 }
						 
						 Term2 = Math.min(qmax.get(j), Nha[ii][j]*vf.get(j));
						 
						 // Explicit Solution
						 Double s_soln = 0.0; Double r_soln = 0.0;
						 if (imputeOR.get(j) && imputeFR.get(j) || Override){
							 
							 Double dummyMax = Math.max(0, Term2 - FlowBet[ii][j]);
							 dummyMax = Math.max(dummyMax, Term2 - c[ii][j]);
							 dummyMax = Math.min(dummyMax, Term2);
							 s_soln = Math.min(Term2, dummyMax);
							 r_soln = c[ii][j] - Term2 + s_soln;
							 
						 } else if(!imputeOR.get(j)){
							 
							 Double s_bnd = Math.max(0, Term2 - c[ii][j]);
							 MultFac = Math.abs(FlowBet[ii][j]/Math.max(100*simulationTimeStep, orFlow_Giv[ii][j]))/2;
							 Double s_vertex = 0.0;
							 if (MultFac<1){
								 s_vertex = Term2 - FlowBet[ii][j];
							 } else {
								 s_vertex = Term2 + orFlow_Giv[ii][j] - c[ii][j];
							 }
							 s_soln = Math.min(Math.max(s_vertex, s_bnd), Term2);
							 r_soln = c[ii][j] - Term2 + s_soln;
							 
						 } else {
							 
							 Double s_bnd = Math.max(0, Term2 - c[ii][j]);
							 MultFac = Math.abs(FlowBet[ii][j]/Math.max(100*simulationTimeStep, frFlow_Giv[ii][j]))/2;
							 Double s_vertex = 0.0;
							 if (MultFac<1){
								 s_vertex = Term2 - FlowBet[ii][j];
							 } else {
								 s_vertex = frFlow_Giv[ii][j];
							 }
							 s_soln = Math.min(Math.max(s_vertex, s_bnd), Term2);
							 r_soln = c[ii][j] - Term2 + s_soln;	
							 
						 }
						 
						 dj1[ii][j] = r_soln;
						 BETA1[ii][j] = s_soln/Term2;
						 OrInp[ii][j] = r_soln;
						 dprev[j] = 0.0;		 
						 
					 }
					 					 
				}
			}
		} // line 758
		
		Double[][] BETAF = new Double[BETA.length][BETA[0].length];
		for (int row = 0;row<BETAF.length;row++){
			for (int col = 0;col<BETAF[0].length;col++){
				BETAF[row][col] = Math.max(0.0, Math.min(BETA[row][col]+BETA1[row][col], 1.0));
			}
		}
		Double[][] DJ = new Double[dj1.length][dj1[0].length];
		for (int row = 0;row<DJ.length;row++){
			for (int col = 0;col<DJ[0].length;col++){
				DJ[row][col] = dj1[row][col] + Demand[row][col];
			}
		}
		
		Double[] DJBound = new Double[BETA.length];
		if (downBoundaryCongested) {
			DJBound = MyUtilities.fetchColumn(c, c.length);
			for (int k=0;i<DJBound.length;k++){
				DJBound[k] -= Math.min(qmax.get(cellData.size()-1), Nha[k][cellData.size()-1]*vf.get(cellData.size()-1));
				if (c[k][cellData.size()-1] <= 0.00001 + qmax.get(cellData.size()-1)){
					DJBound[k] = 0.0;
				}
			}
		}
		
		// *********************************************************************************************
		// ******* Using the DJ, run the model to get onramp flows *************************************
		// *********************************************************************************************
		
		Double[][] Nh_sim = MyUtilities.zerosMatrix(STime.length, cellData.size());
		
		Double[][] qj = MyUtilities.zerosMatrix(STime.length, cellData.size());
		Double[] InFlow = MyUtilities.scaleVector(MyUtilities.onesVector(cellData.size()),0.0);
		Double[] OutFlow = MyUtilities.scaleVector(MyUtilities.onesVector(cellData.size()),0.0);
		
		Double InQ = 0.0;
		MyUtilities.assignRow(Nh_sim, MyUtilities.fetchRow(measuredDensity,1), 0);
		Arrays.fill(dprev,0.0);
		Double[][] OnrampInput = MyUtilities.zerosMatrix(STime.length, cellData.size());
		Double BoundDprev = 0.0;
		Double[] dummyVect = new Double[cellData.size()]; 
		Double[] BoundaryRampFlow = new Double[STime.length];
		Double[] BoundaryRampInp = new Double[STime.length];
			
		for (int ii=1;ii<STime.length;ii++){ // line 786 
			
			dummyVect = MyUtilities.fetchRow(Nh_sim, ii-1);
			dummyVect = MyUtilities.addVectors(dummyVect, InFlow);
			dummyVect = MyUtilities.addVectors(dummyVect, MyUtilities.scaleVector(OutFlow, -1.0));
			Nh_sim = MyUtilities.assignRow(Nh_sim, dummyVect, ii);
			
			OutFlow[OutFlow.length-1] = Math.min(Nh_sim[ii][Nh_sim[0].length-1]*vf.get(vf.size()-1), qmax.get(qmax.size()-1));
			
			if (downBoundaryCongested){
				
				if (DJBound[ii]+OutFlow[OutFlow.length-1]>qmax.get(qmax.size()-1)){
					
					BoundaryRampFlow[ii] = qmax.get(qmax.size()-1)/(DJBound[ii]+OutFlow[OutFlow.length-1])*DJBound[ii];
					BoundaryRampInp[ii] = Math.max(0, DJBound[ii] - BoundDprev);
					
					BoundDprev = DJBound[ii] - BoundaryRampFlow[ii];
					OutFlow[OutFlow.length-1] = qmax.get(qmax.size()-1)/(DJBound[ii]+OutFlow[OutFlow.length-1])*OutFlow[OutFlow.length-1];
					
				} else {
					
					BoundaryRampFlow[ii] = DJBound[ii];
					BoundaryRampInp[ii] = Math.max(0, DJBound[ii] - BoundDprev);
					BoundDprev = DJBound[ii] - BoundaryRampFlow[ii];
					
				}
				
			}
				
				InQ += inputFLW[ii];
				
				Double Capacity = w.get(0)*(rhojam.get(0)-Nh_sim[ii][0]);
				InFlow[0] = Math.min(Capacity, InQ);
				InQ -= InFlow[0];
				
				for (int j=1;j<numberOfNodes-1;j++){
					
					Double[][] Beta = {{(double) Math.round((1-BETAF[ii][j-1])*10000)/10000,1.0},{(double) Math.round(BETAF[ii][j-1]*10000)/10000,0.0}};
					Double[] Demands = {Math.min(qmax.get(j-1), Nh_sim[ii][j-1]*vf.get(j-1)), DJ[ii][j-1]};
					
					Double[] DemandOut = {Beta[0][0]*Demands[0]+Beta[0][1]*Demands[1],Beta[1][0]*Demands[0]+Beta[1][1]*Demands[1]}; // 2x2 times 2x1 matrix multiplication
					Double[] Capacities = {Math.min(qmax.get(j), w.get(j)*(rhojam.get(j)-Nh_sim[ii][j])),10000.0};
					
					Double[][] Dij = new Double[Beta.length][Beta[0].length]; // line 820
					Dij[0][0] = Beta[0][0]*Demands[0]; Dij[0][1] = Beta[0][1]*Demands[1];
					Dij[1][0] = Beta[1][0]*Demands[0]; Dij[1][1] = Beta[1][1]*Demands[1];
					
					if (DemandOut[0] != 0){
						Dij[0][0] *= Math.min(DemandOut[0], Capacities[0])/DemandOut[0];
						Dij[1][0] *= Math.min(DemandOut[0], Capacities[0])/DemandOut[0];
					}
					
					if (DemandOut[1] != 0){
						Dij[0][1] *= Math.min(DemandOut[1], Capacities[1])/DemandOut[1];
						Dij[1][1] *= Math.min(DemandOut[1], Capacities[1])/DemandOut[1];
					}
					
					Demands[0] = Dij[0][0]+Dij[1][0]; Demands[1] = Dij[0][1]+Dij[1][1];
					
					Double[] AdjFact = new Double[2];
					Double term1 = Dij[0][0]/(Demands[0]*Beta[0][0]);
					Double term2 = Dij[1][0]/(Demands[0]*Beta[1][0]);
					Double term3 = Dij[0][1]/(Demands[1]*Beta[0][1]);
					Double term4 = Dij[1][1]/(Demands[1]*Beta[1][1]);
					
					if (!term1.equals(java.lang.Double.NaN) && !term2.equals(java.lang.Double.NaN)){
						AdjFact[0] = Math.min(term1, term2);
					} else if (term1.equals(java.lang.Double.NaN)){
						AdjFact[0] = term2;
					} else {
						AdjFact[0] = term1;
					}
					
					if (!term3.equals(java.lang.Double.NaN) && !term4.equals(java.lang.Double.NaN)){
						AdjFact[1] = Math.min(term3, term4);
					} else if (term3.equals(java.lang.Double.NaN)){
						AdjFact[1] = term4;
					} else {
						AdjFact[1] = term3;
					}
										
					if (AdjFact[0].equals(java.lang.Double.NaN)){ AdjFact[0] = 0.0; }
					if (AdjFact[1].equals(java.lang.Double.NaN)){ AdjFact[1] = 0.0; }
					
					Double[] flow = new Double[2]; flow[0] = Demands[0]*AdjFact[0]; flow[1] = Demands[1]*AdjFact[1];
					
					OutFlow[j-1] = flow[0];
					OrFlow[ii][j-1] = flow[1];
					
					qj[ii][j-1] = qj[ii][j-1]-flow[1];
					InFlow[j] = Beta[0][0]*flow[0] + Beta[0][1]*flow[1];
					FrFlow[ii][j-1] = Beta[1][0]*flow[0] + Beta[1][1]*flow[1];
					FlowBet[ii][j-1] = InFlow[j]-flow[1];
					OnrampInput[ii][j-1] = DJ[ii][j-1]-dprev[j-1];
					dprev[j-1] = DJ[ii][j-1]-flow[1];
					
				}
				
			} // line 842

		Double[][] dummy1 = MyUtilities.addMatrices(measuredDensity, MyUtilities.scaleMatrix(Nh_sim, -1.0));
		dummy1 = MyUtilities.matrixAbsValue(dummy1);
		Double[] dummy2 = MyUtilities.meanColumns(dummy1);
		Double[] dummy3 = MyUtilities.meanColumns(measuredDensity);
		Double dummy4 = MyUtilities.meanVector(dummy2)/MyUtilities.meanVector(dummy3);
		
		System.out.println(" Final Density Error: " + dummy4*100);
		
		// *********************************************************************************************
		// ******* Check Onramp flows calculated in the previous section *******************************
		// *********************************************************************************************
		
		Double[][] Nh_check = MyUtilities.zerosMatrix(STime.length, cellData.size());
		qj = MyUtilities.zerosMatrix(STime.length, cellData.size());
		InFlow = MyUtilities.scaleVector(MyUtilities.onesVector(cellData.size()),0.0);
		OutFlow = MyUtilities.scaleVector(MyUtilities.onesVector(cellData.size()),0.0);
		Double[][] ORINP = OnrampInput;
		for (int j=0;j<ORINP.length;j++){
			for (int k=0;k<ORINP[0].length;k++){
				ORINP[j][k] = Math.max(0, ORINP[j][k]);
			}
		}
		
		Double[] DJBoundSav = new Double[DJBound.length];
		Double BoundaryQ = 0.0;
		InQ = 0.0;
		Nh_check = MyUtilities.assignRow(Nh_check, measuredDensity[0], 0);
		
		Double[][] InFl = MyUtilities.zerosMatrix(STime.length, cellData.size());
		Double[][] OutFl = MyUtilities.zerosMatrix(STime.length, cellData.size());
	
		for(int ii=1;ii<STime.length;ii++){
			
			InFl = MyUtilities.assignRow(InFl, InFlow, ii);
			OutFl = MyUtilities.assignRow(OutFl, OutFlow, ii);
			Double[] dummyV = MyUtilities.addVectors(MyUtilities.fetchRow(Nh_check, ii-1), InFlow);
			dummyV = MyUtilities.addVectors(dummyV, MyUtilities.scaleVector(OutFlow, -1.0));
			Nh_check = MyUtilities.assignRow(Nh_check, dummyV, ii);
			
			OutFlow[OutFlow.length-1] = Math.min(Nh_check[ii][Nh_check[0].length-1]*vf.get(vf.size()-1), qmax.get(qmax.size()-1));
			
			if(downBoundaryCongested){
				
				BoundaryQ = BoundaryQ + BoundaryRampInp[ii];
				DJBoundSav[ii] = BoundaryQ;
				
				if (downBoundaryCongested){
					
					BoundaryRampFlow[ii] = qmax.get(qmax.size()-1)/(BoundaryQ+OutFlow[OutFlow.length-1])*BoundaryQ;
					OutFlow[OutFlow.length-1] = qmax.get(qmax.size()-1)/(BoundaryQ+OutFlow[OutFlow.length-1])*OutFlow[OutFlow.length-1];
					BoundaryQ = BoundaryQ - BoundaryRampFlow[ii];
					
				} else {
					
					BoundaryRampFlow[ii] = BoundaryQ;
					BoundaryQ = BoundaryQ - BoundaryRampFlow[ii];
					
				}
				
			}
			
			qj = MyUtilities.assignRow(qj, MyUtilities.addVectors(MyUtilities.fetchRow(qj, ii-1), MyUtilities.fetchRow(ORINP, ii)), ii);
			
			InQ = InQ + inputFLW[ii];
			Double Capacity = w.get(0)*(rhojam.get(0)-Nh_check[ii][1]);
			InFlow[0] = Math.min(Capacity, InQ);
			InQ = InQ - InFlow[0];
			
			for (int j=1;j<numberOfNodes-1;j++){
				
				Double[][] Beta = {{(double) Math.round((1-BETAF[ii][j-1])*10000)/10000,1.0},{(double) Math.round(BETAF[ii][j-1]*10000)/10000,0.0}};
				Double[] Demands = {Math.min(qmax.get(j-1), Nh_check[ii][j-1]*vf.get(j-1)), qj[ii][j-1]};
				
				Double[] DemandOut = {Beta[0][0]*Demands[0]+Beta[0][1]*Demands[1],Beta[1][0]*Demands[0]+Beta[1][1]*Demands[1]}; // 2x2 times 2x1 matrix multiplication
				Double[] Capacities = {Math.min(qmax.get(j), w.get(j)*(rhojam.get(j)-Nh_check[ii][j])),10000.0};
				
				Double[][] Dij = new Double[Beta.length][Beta[0].length];
				Dij[0][0] = Beta[0][0]*Demands[0]; Dij[0][1] = Beta[0][1]*Demands[1];
				Dij[1][0] = Beta[1][0]*Demands[0]; Dij[1][1] = Beta[1][1]*Demands[1];
				
				if (DemandOut[0] != 0){
					Dij[0][0] *= Math.min(DemandOut[0], Capacities[0])/DemandOut[0];
					Dij[1][0] *= Math.min(DemandOut[0], Capacities[0])/DemandOut[0];
				}
				
				if (DemandOut[1] != 0){
					Dij[0][1] *= Math.min(DemandOut[1], Capacities[1])/DemandOut[1];
					Dij[1][1] *= Math.min(DemandOut[1], Capacities[1])/DemandOut[1];
				}
				
				Demands[0] = Dij[0][0]+Dij[1][0]; Demands[1] = Dij[0][1]+Dij[1][1];

				Double[] AdjFact = new Double[2];
				Double term1 = Dij[0][0]/(Demands[0]*Beta[0][0]);
				Double term2 = Dij[1][0]/(Demands[0]*Beta[1][0]);
				Double term3 = Dij[0][1]/(Demands[1]*Beta[0][1]);
				Double term4 = Dij[1][1]/(Demands[1]*Beta[1][1]);
				
				if (!term1.equals(java.lang.Double.NaN) && !term2.equals(java.lang.Double.NaN)){
					AdjFact[0] = Math.min(term1, term2);
				} else if (term1.equals(java.lang.Double.NaN)){
					AdjFact[0] = term2;
				} else {
					AdjFact[0] = term1;
				}
				
				if (!term3.equals(java.lang.Double.NaN) && !term4.equals(java.lang.Double.NaN)){
					AdjFact[1] = Math.min(term3, term4);
				} else if (term3.equals(java.lang.Double.NaN)){
					AdjFact[1] = term4;
				} else {
					AdjFact[1] = term3;
				}
									
				if (AdjFact[0].equals(java.lang.Double.NaN)){ AdjFact[0] = 0.0; }
				if (AdjFact[1].equals(java.lang.Double.NaN)){ AdjFact[1] = 0.0; }
				
				Double[] flow = new Double[2]; flow[0] = Demands[0]*AdjFact[0]; flow[1] = Demands[1]*AdjFact[1];
				
				
				OutFlow[j-1] = flow[0];
				OrFlow[ii][j-1] = flow[1];
				
				qj[ii][j-1] = qj[ii][j-1]-flow[1];
				InFlow[j] = Beta[0][0]*flow[0] + Beta[0][1]*flow[1];
				FrFlow[ii][j-1] = Beta[1][0]*flow[0] + Beta[1][1]*flow[1];
				FlowBet[ii][j-1] = InFlow[j]-flow[1];
								
			}
			
		}
		
		// calculate density and flow errors at the end of the method:
		
		// Final Density Error (Matlab line: DensityError=num2str(mean(mean(abs(Density-Nh)))./mean(mean(Nh))*100);)
		Double[][] ddummy1 = MyUtilities.addMatrices(measuredDensity, MyUtilities.scaleMatrix(Nh_sim, -1.0));
		ddummy1 = MyUtilities.matrixAbsValue(ddummy1);
		Double[] ddummy2 = MyUtilities.meanColumns(ddummy1);
		Double[] ddummy3 = MyUtilities.meanColumns(measuredDensity);
		Double ddummy4 = MyUtilities.meanVector(ddummy2)/MyUtilities.meanVector(ddummy3);
		
		System.out.println(" Final Density Error: " + ddummy4*100);
		
		// Final Flow Error (Matlab line: FlowError=num2str(mean(mean(abs(FlowCompare-Flow)))./mean(mean(Flow))*100);)
		Double[][] FlowCompare = MyUtilities.zerosMatrix(Nh_check.length, Nh_check[0].length);
		FlowCompare = MyUtilities.assignColumn(FlowCompare, MyUtilities.fetchColumn(InFl, 0), 0);
		for (int j=1;j<FlowCompare[0].length;j++){
			MyUtilities.assignColumn(FlowCompare, MyUtilities.fetchColumn(FlowBet, j-1), j);
		}
		Double[][] dddummy1 = MyUtilities.addMatrices(measuredFlow, MyUtilities.scaleMatrix(FlowCompare, -1.0));
		dddummy1 = MyUtilities.matrixAbsValue(dddummy1);
		Double[] dddummy2 = MyUtilities.meanColumns(dddummy1);
		Double[] dddummy3 = MyUtilities.meanColumns(measuredFlow);
		Double dddummy4 = MyUtilities.meanVector(dddummy2)/MyUtilities.meanVector(dddummy3);
		
		System.out.println(" Final Flow Error: " + dddummy4*100);
		
		int tratio = (int) (demandTimeStep / simulationTimeStep);
		// Fill in the missing fields in the cells:
//		for (int j=0;j<this.cellData.size();j++){
//			Double[] hand = MyUtilities.meanRows(MyUtilities.reshapeVectorIntoMatrix(MyUtilities.fetchColumn(OnrampInput,j),tratio,(int) (24/this.demandTimeStep)));
//			ArrayList<Double> cardsList = new ArrayList<Double>(Arrays.asList(hand));
//			cellData.get(j+1).setOnRampInput();
//		}
//		
			
	} // end of method run()
	
	private void initializeDataMatrices(){
		
		measuredDensity = new Double[detectorList.values().iterator().next().getDensityData().size()][cellData.size()];
		measuredFlow = new Double[detectorList.values().iterator().next().getDensityData().size()][cellData.size()];
		measuredSpeed = new Double[detectorList.values().iterator().next().getDensityData().size()][cellData.size()];
		measuredOrFlow = new Double[detectorList.values().iterator().next().getDensityData().size()][cellData.size()];
		measuredFrFlow = new Double[detectorList.values().iterator().next().getDensityData().size()][cellData.size()];
		
	}

}
