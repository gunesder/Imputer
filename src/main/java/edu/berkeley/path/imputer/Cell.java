package edu.berkeley.path.imputer;

import java.util.ArrayList;
import java.util.LinkedList;

import org.apache.commons.math.linear.BlockRealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.joda.time.Interval;

public class Cell {
	
	// fields
	private LinkedList<Link> links = new LinkedList<Link>();
	private Detector detectorML;
	private Detector detectorHOV;
	private boolean majorOR;
	private boolean majorFR;
	private ArrayList<Double> orFlowNonImputed = new ArrayList<Double>();
	private ArrayList<Double> frFlowNonImputed = new ArrayList<Double>();
	private ArrayList<Boolean> imputeOR = new ArrayList<Boolean>();
	private ArrayList<Boolean> imputeFR = new ArrayList<Boolean>();
	private ArrayList<Double> outFLow = new ArrayList<Double>();
	private ArrayList<Double> flowCompare = new ArrayList<Double>();
	private ArrayList<Double> simDensity = new ArrayList<Double>();
	private ArrayList<Double> onRampFlow = new ArrayList<Double>();
	private ArrayList<Double> offRampFlow = new ArrayList<Double>();
	private ArrayList<Double> beta = new ArrayList<Double>();
	private ArrayList<Double> onRampInput = new ArrayList<Double>();
	private BlockRealMatrix measuredOnrampFlow = new BlockRealMatrix(288,1); // static pre-allocation, might fix later
	private BlockRealMatrix measuredOfframpFlow = new BlockRealMatrix(288,1); // static pre-allocation, might fix later
	private ArrayList<Double> Velocity = new ArrayList<Double>();
	private ArrayList<Double> Demand = new ArrayList<Double>();
	private ArrayList<Integer> onrampsPerLink = new ArrayList<Integer>();
	private ArrayList<Integer> offrampsPerLink = new ArrayList<Integer>();
	
	// getters and setters
	public LinkedList<Link> getLinks() {
		return links;
	}
	public BlockRealMatrix getMeasuredOnrampFlow() {
		return measuredOnrampFlow;
	}
	public void setMeasuredOnrampFlow(BlockRealMatrix measuredOnrampFlow) {
		this.measuredOnrampFlow = measuredOnrampFlow;
	}
	public BlockRealMatrix getMeasuredOfframpFlow() {
		return measuredOfframpFlow;
	}
	public void setMeasuredOfframpFlow(BlockRealMatrix measuredOfframpFlow) {
		this.measuredOfframpFlow = measuredOfframpFlow;
	}
	public ArrayList<Integer> getOnrampsPerLink() {
		return onrampsPerLink;
	}
	public void setOnrampsPerLink(ArrayList<Integer> onrampsPerLink) {
		this.onrampsPerLink = onrampsPerLink;
	}
	public ArrayList<Integer> getOfframpsPerLink() {
		return offrampsPerLink;
	}
	public void setOfframpsPerLink(ArrayList<Integer> offrampsPerLink) {
		this.offrampsPerLink = offrampsPerLink;
	}
	public void setLinks(LinkedList<Link> links) {
		this.links = links;
	}
	public Detector getDetectorML() {
		return detectorML;
	}
	public void setDetectorML(Detector detectorML) {
		this.detectorML = detectorML;
	}
	public Detector getDetectorHOV() {
		return detectorHOV;
	}
	public void setDetectorHOV(Detector detectorHOV) {
		this.detectorHOV = detectorHOV;
	}
	public boolean isMajorOR() {
		return majorOR;
	}
	public void setMajorOR(boolean majorOR) {
		this.majorOR = majorOR;
	}
	public boolean isMajorFR() {
		return majorFR;
	}
	public void setMajorFR(boolean majorFR) {
		this.majorFR = majorFR;
	}
	public ArrayList<Double> getOrFlowNonImputed() {
		return orFlowNonImputed;
	}
	public void setOrFlowNonImputed(ArrayList<Double> orFlowNonImputed) {
		this.orFlowNonImputed = orFlowNonImputed;
	}
	public ArrayList<Double> getFrFlowNonImputed() {
		return frFlowNonImputed;
	}
	public void setFrFlowNonImputed(ArrayList<Double> frFlowNonImputed) {
		this.frFlowNonImputed = frFlowNonImputed;
	}
	public ArrayList<Boolean> isImputeOR() {
		return imputeOR;
	}
	public void setImputeOR(ArrayList<Boolean> imputeOR) {
		this.imputeOR = imputeOR;
	}
	public ArrayList<Boolean> isImputeFR() {
		return imputeFR;
	}
	public void setImputeFR(ArrayList<Boolean> imputeFR) {
		this.imputeFR = imputeFR;
	}
	public ArrayList<Double> getOutFLow() {
		return outFLow;
	}
	public void setOutFLow(ArrayList<Double> outFLow) {
		this.outFLow = outFLow;
	}
	public ArrayList<Double> getFlowCompare() {
		return flowCompare;
	}
	public void setFlowCompare(ArrayList<Double> flowCompare) {
		this.flowCompare = flowCompare;
	}
	public ArrayList<Double> getSimDensity() {
		return simDensity;
	}
	public void setSimDensity(ArrayList<Double> simDensity) {
		this.simDensity = simDensity;
	}
	public ArrayList<Double> getOnRampFlow() {
		return onRampFlow;
	}
	public void setOnRampFlow(ArrayList<Double> onRampFlow) {
		this.onRampFlow = onRampFlow;
	}
	public ArrayList<Double> getOffRampFlow() {
		return offRampFlow;
	}
	public void setOffRampFlow(ArrayList<Double> offRampFlow) {
		this.offRampFlow = offRampFlow;
	}
	public ArrayList<Double> getBeta() {
		return beta;
	}
	public void setBeta(ArrayList<Double> beta) {
		this.beta = beta;
	}
	public ArrayList<Double> getOnRampInput() {
		return onRampInput;
	}
	public void setOnRampInput(ArrayList<Double> onRampInput) {
		this.onRampInput = onRampInput;
	}	
	public BlockRealMatrix getmeasuredOnrampFlow() {
		return measuredOnrampFlow;
	}
	public void setmeasuredOnrampFlow(BlockRealMatrix measuredOnrampFlow) {
		this.measuredOnrampFlow = measuredOnrampFlow;
	}
	public BlockRealMatrix getmeasuredOfframpFlow() {
		return measuredOfframpFlow;
	}
	public void setmeasuredOfframpFlow(BlockRealMatrix measuredOfframpFlow) {
		this.measuredOfframpFlow = measuredOfframpFlow;
	}
	public ArrayList<Double> getVelocity() {
		return Velocity;
	}
	public void setVelocity(ArrayList<Double> velocity) {
		Velocity = velocity;
	}
	public ArrayList<Double> getDemand() {
		return Demand;
	}
	public void setDemand(ArrayList<Double> demand) {
		Demand = demand;
	}
	
	// constructors
	
	public Cell(int datasize){
		measuredOnrampFlow = new BlockRealMatrix(datasize,1);
		measuredOfframpFlow = new BlockRealMatrix(datasize,1);
	}
	
	// methods
	public void addLink(Link x){
		this.links.add(x);
	}
	
	public void addToOnrampPerLink(int i){
		this.onrampsPerLink.add(i);
	}
	
	public void addToOfframpPerLink(int i){
		this.offrampsPerLink.add(i);
	}
	
	public void addToImputeOR(boolean b){
		this.imputeOR.add(b);
	}
	
	public void addToImputeFR(boolean b){
		this.imputeFR.add(b);
	}
	
	public void addColumnToMeasuredOnrampFlow (RealVector v){
		if (this.measuredOnrampFlow.getColumnDimension() == 1){
			this.measuredOnrampFlow.setColumnVector(0, v);
		} else {
			int newColumnDimension = this.measuredOnrampFlow.getColumnDimension() + 1;
			BlockRealMatrix dummy = new BlockRealMatrix(links.getFirst().getDetectorML().getDensityData().size(),newColumnDimension);
			for (int i = 0;i<newColumnDimension - 1;i++){
				dummy.setColumn(i,this.measuredOnrampFlow.getColumn(i));
			}
			dummy.setColumnVector(dummy.getColumnDimension(), v);
			this.measuredOnrampFlow = dummy;
		}
	}
	
	public void addColumnToMeasuredOfframpFlow (RealVector v){
		if (this.measuredOnrampFlow.getColumnDimension() == 1){
			this.measuredOnrampFlow.setColumnVector(0, v);
		} else {
			int newColumnDimension = this.measuredOfframpFlow.getColumnDimension() + 1;
			BlockRealMatrix dummy = new BlockRealMatrix(links.getFirst().getDetectorML().getDensityData().size(),newColumnDimension);
			for (int i = 0;i<newColumnDimension - 1;i++){
				dummy.setColumn(i,this.measuredOfframpFlow.getColumn(i));
			}
			dummy.setColumnVector(dummy.getColumnDimension(), v);
			this.measuredOfframpFlow = dummy;
		}
	}
	
}
