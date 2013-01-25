package edu.berkeley.path.imputer;

import java.util.*;

import lpsolve.LpSolve;
import lpsolve.LpSolveException;

public class ImputationCoreAlgorithm {
	
	// fields
		// primary fields
		private LinkedList<Cell> cellData = new LinkedList<Cell>();
		private HashMap<Integer,Detector> detectorList = new HashMap<Integer,Detector>();
		private double demandTimeStep = 5.0/60.0; // [hours] 5 minute default
		private double simulationTimeStep = 5.0/60.0/60.0; // [hours] 5 second default
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
		private double[][] measuredDensity = new double[289][1];
		private double[][] measuredFlow = new double[289][1];
		private double[][] measuredSpeed = new double[289][1];
		private double[][] measuredOrFlow = new double[289][1];
		private double[][] measuredFrFlow = new double[289][1];		
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
		for (int j=0;j<imputeOR.size();j++){
			impute.add(imputeOR.get(j)|imputeFR.get(j));
		}
		
		int numberOfNodes = measuredDensity[0].length + 1;
		
		double[][] BETA = new double[FrFlow.length][FrFlow[0].length];
		double[][] dj = new double[FrFlow.length][FrFlow[0].length];
		double[] zeroVector = new double[FrFlow.length];
		Arrays.fill(zeroVector, 0);
		BETA = MyUtilities.appendColumn(BETA, zeroVector);
		dj = MyUtilities.appendColumn(dj, zeroVector);
		double[][] OrINP = new double[dj.length][dj[0].length];
		double[] dprev = new double[numberOfNodes-2];
		Arrays.fill(dprev, 0);
		double[][] ErrA = new double[BETA.length][BETA[0].length];
		double[][] ErrB = new double[BETA.length][BETA[0].length];
		
		// skipped line 198 for now
		
		// *********************************************************************************************
		// ************** Derive dj and BETA from known ramp flows *************************************
		// *********************************************************************************************
		
		double[] STime = MyUtilities.createIncrementVector(0, 24-this.simulationTimeStep, this.simulationTimeStep);
		for (int ii=1;ii<STime.length;ii++){
			
			for (int j=0;j<numberOfNodes-2;j++){
				
				if (!impute.get(j)){
					
					 double Term1 = this.qmax.get(j+1) < this.w.get(j+1)*(this.rhojam.get(j+1) - measuredDensity[ii][j+1]) ? this.qmax.get(j+1) : this.w.get(j+1)*(this.rhojam.get(j+1)-measuredDensity[ii][j+1]);
					 double Term2 = this.qmax.get(j) < measuredDensity[ii][j]*this.vf.get(j) ? this.qmax.get(j) : measuredDensity[ii][j]*this.vf.get(j);
					 
					 double r = OrFlow[ii][j];
					 double s = FrFlow[ii][j] < Term2 ? FrFlow[ii][j] : Term2;
					 
					 if (Term2+r-s>Term1){
						 
						 double a1 = Term1 - r;
						 double a2 = Term2*r;
						 double c1 = a2;
						 double b1 = -s;
						 double b2 = Term2*(Term1+s);
						 double c2 = s*Term2;
						 int[] colno = new int[4];
				         double[] row = new double[4];
						 					 
						 try {
							 
							LpSolve solver = LpSolve.makeLp(0, 4); // 0 constraints, 4 variables to start with
							solver.strSetObjFn("1 1 0 0"); // cost is the sum of the first two variables
							
							solver.setAddRowmode(true);  /* makes building the model faster if it is done rows by row */

				            /* construct first row (-x1 + a1x3 +a2x4 <= c1) */
				            j = 0;

				            colno[j] = 1; /* first column */
				            row[j++] = -1;

				            colno[j] = 2; /* second column */
				            row[j++] = 0;
				            
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
						dprev[j] = 0;
						OrINP[ii][j] = dj[ii][j] - dprev[j];
						BETA[ii][j] = s/(0.0000000001 > Term2 ? 0.0000000001 : Term2);
					 }
					 
				}
			}
		}
		
		double[][] Demand = dj; 
		
		// *********************************************************************************************
		// ************** Check that ramps are present whenever imputation is enabled ******************
		// *********************************************************************************************
		
		double[] lowerBounds = new double[numberOfNodes-2];
		double[] upperBounds = new double[numberOfNodes-2];
		for (int ind = 0;ind < numberOfNodes - 2; ind++){
			if (frPresent.get(ind)) lowerBounds[ind] = 0; else lowerBounds[ind] = 1;
			if (orPresent.get(ind)) upperBounds[ind] = 0; else upperBounds[ind] = 1;
			if (!frPresent.get(ind) & !orPresent.get(ind) & impute.get(ind)){
				
			}			
		}
		if (downBoundaryCongested){
			lowerBounds[numberOfNodes-2] = 1;
			upperBounds[numberOfNodes-2] = 0;
		}
		
		// *********************************************************************************************
		// ************************* Learning Algorithm ************************************************
		// *********************************************************************************************
		
		/* Learning algorithm front matter, gains, user-defined settings, etc. (hard coded for now) */
		double GM = 40; // these are user-defined gains for the adaptive learning algorithm and are not necessarily all very sensitive
		double G1 = 1*GM;
		double G2 = 0.0001*GM;
		double percTol = 0.001;
		double percTol2 = 0.2;
		int maxLim = 100;
		int iterMax = 25;
		int[] iterTrigger = {5,9,12,16,20}; // iteration indices when the trigger algorithm kicks in
		double derivativeBound = 2;
		double startBound = 1;
		
		// initialize c: this is the effective demand vector into each cell
		double[][] c = new double[STime.length][cellData.size()-1];
		double[] dummy = new double[measuredDensity.length];
		for (int j=1;j<numberOfNodes-1;j++){
			dummy = MyUtilities.scaleVector(MyUtilities.onesVector(measuredDensity.length), rhojam.get(j));
			dummy = MyUtilities.addVectors(dummy,MyUtilities.scaleVector(MyUtilities.fetchColumn(measuredDensity, j),-1));
			c = MyUtilities.assignColumn(c, MyUtilities.scaleVector(dummy, 1.5*w.get(j)), j);
		}
		
		if (downBoundaryCongested){
			MyUtilities.assignColumn(c, MyUtilities.scaleVector(MyUtilities.onesVector(c.length), qmax.get(numberOfNodes-1)), numberOfNodes-1);
			impute.set(numberOfNodes, true);
		}
		
		double[][] cBest = c;
		int iterBound_D_Beta = iterMax+1;
		int startIterBound = 10;
		
		boolean[] flags = new boolean[numberOfNodes];
		
		// Learning Algorithm Loop
		for (int iter = 1;iter <= iterMax;iter++){
			
			flags[0] = true;
			
			double[][] mode = MyUtilities.zerosMatrix(c.length, c[0].length-1);
			
			double InQ = 0;
			double[][] Nh = MyUtilities.zerosMatrix(STime.length, cellData.size());
			Nh = MyUtilities.assignRow(Nh, measuredDensity[0], 0);
			double[] cj = c[0];
			
			for (int k=0;k<STime.length-1;k++){ // starts at line 328 and ends at line 584
				
				double[] Limit = new double[cellData.size()];
				for (int j=0;j<Limit.length;j++){
					Limit[j] = qmax.get(j) < w.get(j)*(rhojam.get(j)-Nh[k][j]) ? qmax.get(j) : w.get(j)*(rhojam.get(j)-Nh[k][j]);
				}
				
				double[] cjprev = cj;
				cj = c[k];
				
				/* This whole section does what line 333 does in Matlab, may need to make into a utility method if it is needed somewhere else*/
				double[] dummy1 = MyUtilities.subtractVectors(cjprev,MyUtilities.scaleVector(MyUtilities.onesVector(cjprev.length), derivativeBound));
				double[] dummy2 = MyUtilities.addVectors(cjprev, MyUtilities.scaleVector(MyUtilities.onesVector(cjprev.length), derivativeBound));
				double[] dummy3 = new double[cj.length];
				double[] dummy4 = new double[cj.length];
				
				for (i=0;i<cj.length;i++){
					dummy3[i] = cj[i] < dummy2[i] ? cj[i] : dummy2[i];
				}
				
				for (i=0;i<cj.length;i++){
					dummy4[i] = dummy1[i] > dummy3[i] ? dummy1[i] : dummy3[i];
				}
				/* This whole section does what line 333 does in Matlab, may need to make into a utility method if it is needed somewhere else*/
				
				for (int j=0;j<numberOfNodes-2;j++){
					
					double NjVfj = qmax.get(j) < Nh[k][j]*vf.get(j) ? qmax.get(j) : Nh[k][j]*vf.get(j);
					
					if(!impute.get(j)){
						cj[j] = NjVfj*(1-BETA[k][j]) + Demand[k][j];
					} else {
						if (iter>startIterBound){
							if (lowerBounds[j] !=0){
								cj[j] = cj[j] > NjVfj ? cj[j] : NjVfj;
							}
							if (upperBounds[j] != 0){
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
				double inFlow = Limit[0] < InQ ? Limit[0] : InQ;
				InQ = InQ - inFlow;
				
				boolean needBoundaryImpute = downBoundaryCongested & (boundaryVelocity[k]<boundaryVF*0.9);
				if (needBoundaryImpute){
					flags[numberOfNodes-1] = true;
				}
				
				int j=0;
				while (j<numberOfNodes-1){ // starts at line 374 and ends at line 580
					j = j+1;
					double NjVfj = Nh[k][j]*vf.get(j) < qmax.get(j) ? Nh[k][j]*vf.get(j) : qmax.get(j);
					double Njm1Vfjm1 = 0; // never used, just a workaround. It is assigned a new value for each node except for the first one, which doesn't undergo imputation anyway
					if (j>1){
						Njm1Vfjm1 = qmax.get(j-1) < Nh[k][j-1]*vf.get(j-1) ? qmax.get(j-1) : Nh[k][j-1]*vf.get(j-1);	
					}
					double Wj = qmax.get(j) < w.get(j)*(rhojam.get(j)-Nh[k][j]) ? qmax.get(j) : w.get(j)*(rhojam.get(j)-Nh[k][j]);
					double Wjp1 = 0;
					if (j == numberOfNodes-1){
						Wjp1 = qmax.get(j);
					} else {
						Wjp1 = qmax.get(j+1) < w.get(j+1)*(rhojam.get(j+1)-Nh[k][j+1]) ? qmax.get(j+1) : w.get(j+1)*(rhojam.get(j+1)-Nh[k][j+1]);
					}
					
					if (!flags[j] && !flags[j+1]){ // starts at line 391 and ends at line 577
						mode[k][j] = 1;
						
						double BoundL = lowerBounds[j-1];
						double BoundU = upperBounds[j-1];
						
						double Nh0 = measuredDensity[k+1][j] - (Nh[k][j]-NjVfj+cj[j-1]);
						double NHt = Nh0/(1+G1);
						
						if (impute.get(j-1)) { // Line 402 to 424
							
							cj[j-1] = cj[j-1] + G1*NHt;
							
							double UBound = 0;
							double LBound = 0;
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
						
					} else if (!flags[j] & flags[j+1]){
						mode[k][j] = 2;
						
						if (impute.get(j) & cj[j] > Wjp1*(1+percTol) & cj[j] < Wjp1*(1+percTol) & (Nh[k][j]+cj[j-1]-NjVfj)>measuredDensity[k+1][j]){
							cj[j] = cj[j] < Wjp1*(1-percTol/100) ? cj[j] : Wjp1*(1-percTol/100);
							flags[j+1] = false;
							j = j-1;
							continue;							
						}
						
						cj[j] = cj[j] > Wjp1 ? cj[j] : Wjp1;
						
						double BoundL1 = lowerBounds[j-1];
						double BoundU1 = upperBounds[j-1];
						double BoundL2 = lowerBounds[j];
						double BoundU2 = upperBounds[j];
						
						double Nh0 = measuredDensity[k+1][j] - (Nh[k][j]+cj[j-1]-Wjp1*NjVfj/cj[j]);
						double NHt = Nh0 / (1+G1*(impute.get(j-1) ? 1 : 0)+G2*Wjp1*NjVfj*(impute.get(j) ? 1 : 0));
						
						if (impute.get(j)){
							cj[j] = 1 / (1/cj[j]-(G2*NHt < 1/cj[j]*0.99 ? G2*NHt : 1/cj[j]*0.99));
							
							double UBound2 = 0;
							if (BoundU2 != 0){
								UBound2 =  (maxLim*Wjp1 < NjVfj ? maxLim*Wjp1 : NjVfj) > 0.01 ? (maxLim*Wjp1 < NjVfj ? maxLim*Wjp1 : NjVfj) : 0.01;
							} else if (iter > iterBound_D_Beta){
								UBound2 = maxLim*Wjp1 < (Wjp1 > NjVfj+maxD.get(j) ? Wjp1 : NjVfj+maxD.get(j)) ? maxLim*Wjp1 : (Wjp1 > NjVfj+maxD.get(j) ? Wjp1 : NjVfj+maxD.get(j));
							} else {
								UBound2 = maxLim*Wjp1;
							}
							
							double LBound2 = NjVfj*BoundL2 > Wjp1 ? NjVfj*BoundL2 : Wjp1;
							
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
							
							double UBound1 = 0;
							double LBound1 = 0;
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
						
					} else if (flags[j] & !flags[j+1]){
						
						mode[k][j] = 3;
						if (j==1){
							Nh[k+1][j] = Nh[k][j] + inFlow - (Nh[k][j]*vf.get(j) < qmax.get(j) ? Nh[k][j]*vf.get(j) : qmax.get(j));
						} else {
							Nh[k+1][j] = Nh[k][j] + (qmax.get(j) < w.get(j)*(rhojam.get(j)-Nh[k][j]) ? qmax.get(j) : w.get(j)*(rhojam.get(j)-Nh[k][j])) - (Nh[k][j]*vf.get(j) < qmax.get(j) ? Nh[k][j]*vf.get(j) : qmax.get(j));
						}
						
					} else if (flags[j] & flags[j+1]){
						mode[k][j] = 4;
						
						if (impute.get(j)){
							double BoundL = lowerBounds[j];
							double BoundU = upperBounds[j];
							
							double Nh0 = 0;
							if (j==1){
								Nh0 = measuredDensity[k+1][j] - (Nh[k][j]+inFlow-Wjp1*NjVfj/cj[j]);
							} else {
								Nh0 = measuredDensity[k+1][j] - (Nh[k][j]+Wj-Wjp1*NjVfj/cj[j]);
							}
							
							if (cj[j]>Wjp1*(1-percTol) & cj[j]<Wjp1*(1+percTol) & Nh0<0){
								cj[j] = cj[j] < Wjp1*(1-percTol/100) ? cj[j] : Wjp1*(1-percTol/100);
								flags[j+1] = false;
								j = j-1;
								continue;							
							}
							
							cj[j] = cj[j] > Wjp1 ? cj[j] : Wjp1;
							
							double NHt = Nh0/(1+G2*Wjp1*NjVfj);
							
							cj[j] = 1/(1/cj[j]-(G2*NHt < 1/cj[j]*0.99 ? G2*NHt : 1/cj[j]*0.99));
							
							double UBound = 0;
							if (BoundU != 0){
								UBound = (maxLim*Wjp1 < NjVfj ? maxLim*Wjp1 : NjVfj) > 0.01 ? (maxLim*Wjp1 < NjVfj ? maxLim*Wjp1 : NjVfj) : 0.01;
							} else if (iter > iterBound_D_Beta){
								UBound = maxLim*Wjp1 < (Wjp1 > NjVfj*maxD.get(j) ? Wjp1 : NjVfj*maxD.get(j)) ? maxLim*Wjp1 : (Wjp1 > NjVfj*maxD.get(j) ? Wjp1 : NjVfj*maxD.get(j));
							} else {
								UBound = maxLim*Wjp1;
							}
							
							double LBound = NjVfj*BoundL > Wjp1 ? NjVfj*BoundL : Wjp1;
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
						
						if (j==1){
							Nh[k+1][j] = Nh[k][j] + inFlow - Wjp1*NjVfj/cj[j];
						} else {
							Nh[k+1][j] = Nh[k][j] + Wj - Wjp1*NjVfj/cj[j];
						}
						
						double NHt = measuredDensity[k+1][j] - Nh[k+1][j];
						
					}
					
				} // end of while loop over nodes
				
				for (i=0;i<cj.length;i++){
					cj[i] = 0.01 > cj[i] ? 0.01 : cj[i];
				}
				
				MyUtilities.assignRow(c, cj, k);
				
			} // end of for loop over nodes
			
			// line 586
			
		} // end of for loop over time
		
		
		
		
		
		
		
			
	} // end of method run()
	
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
