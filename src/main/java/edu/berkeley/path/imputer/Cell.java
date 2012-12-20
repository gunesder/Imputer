package edu.berkeley.path.imputer;

import java.util.ArrayList;
import java.util.LinkedList;

public class Cell {
	
	// fields
	private LinkedList<Link> links = new LinkedList<Link>();
	private Detector detectorML;
	private Detector detectorHOV;
	private boolean majorOR;
	private boolean majorFR;
	private ArrayList<Double> orFlowNonImputed = new ArrayList<Double>();
	private ArrayList<Double> frFlowNonImputed = new ArrayList<Double>();
	private boolean imputeOR;
	private boolean imputeFR;
	private ArrayList<Double> outFLow = new ArrayList<Double>();
	private ArrayList<Double> flowCompare = new ArrayList<Double>();
	private ArrayList<Double> simDensity = new ArrayList<Double>();
	private ArrayList<Double> onRampFlow = new ArrayList<Double>();
	private ArrayList<Double> offRampFlow = new ArrayList<Double>();
	private ArrayList<Double> beta = new ArrayList<Double>();
	private ArrayList<Double> onRampInput = new ArrayList<Double>();
	private ArrayList<Double> givenOnrampFlow = new ArrayList<Double>();
	private ArrayList<Double> givenOfframpFlow = new ArrayList<Double>();
	private ArrayList<Double> Velocity = new ArrayList<Double>();
	private ArrayList<Double> Demand = new ArrayList<Double>();
	
	// getters and setters
	public LinkedList<Link> getLinks() {
		return links;
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
	public boolean isImputeOR() {
		return imputeOR;
	}
	public void setImputeOR(boolean imputeOR) {
		this.imputeOR = imputeOR;
	}
	public boolean isImputeFR() {
		return imputeFR;
	}
	public void setImputeFR(boolean imputeFR) {
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
	public ArrayList<Double> getGivenOnrampFlow() {
		return givenOnrampFlow;
	}
	public void setGivenOnrampFlow(ArrayList<Double> givenOnrampFlow) {
		this.givenOnrampFlow = givenOnrampFlow;
	}
	public ArrayList<Double> getGivenOfframpFlow() {
		return givenOfframpFlow;
	}
	public void setGivenOfframpFlow(ArrayList<Double> givenOfframpFlow) {
		this.givenOfframpFlow = givenOfframpFlow;
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
	
	// methods
	public void addLink(Link x){
		this.links.add(x);
	}
	
}
