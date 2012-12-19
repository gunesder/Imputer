package edu.berkeley.path.imputer;

import java.util.ArrayList;

public class Detector {
	
	// fields
	private int sensorID;
	private CharSequence sensorType;
	private String sourceFileName;
	private String sourceAddress;
	private double healthStatus;
	private ArrayList<Double> speedData = new ArrayList<Double>();
	private ArrayList<Double> densityData = new ArrayList<Double>();
	private ArrayList<Double> flowData = new ArrayList<Double>();
	private int freewayNumber;
	private CharSequence freewayDirection;
	private double absolutePM;
	private double latitude;
	private double longitude;
	private double detectorLength;
	private CharSequence detectorName;
	private int numberOfLanes;
	private int linkAssoc;
	
	// getters and setters
	public int getLinkAssoc() {
		return linkAssoc;
	}
	public void setLinkAssoc(int linkAssoc) {
		this.linkAssoc = linkAssoc;
	}
	public int getFreewayNumber() {
		return freewayNumber;
	}
	public void setFreewayNumber(int freewayNumber) {
		this.freewayNumber = freewayNumber;
	}
	public CharSequence getFreewayDirection() {
		return freewayDirection;
	}
	public void setFreewayDirection(CharSequence charSequence) {
		this.freewayDirection = charSequence;
	}
	public double getAbsolutePM() {
		return absolutePM;
	}
	public void setAbsolutePM(double absolutePM) {
		this.absolutePM = absolutePM;
	}
	public double getLatitude() {
		return latitude;
	}
	public void setLatitude(double latitude) {
		this.latitude = latitude;
	}
	public double getLongitude() {
		return longitude;
	}
	public void setLongitude(double longitude) {
		this.longitude = longitude;
	}
	public double getDetectorLength() {
		return detectorLength;
	}
	public void setDetectorLength(double detectorLength) {
		this.detectorLength = detectorLength;
	}
	public CharSequence getDetectorName() {
		return detectorName;
	}
	public void setDetectorName(CharSequence charSequence) {
		this.detectorName = charSequence;
	}
	public int getNumberOfLanes() {
		return numberOfLanes;
	}
	public void setNumberOfLanes(int numberOfLanes) {
		this.numberOfLanes = numberOfLanes;
	}
	public int getSensorID() {
		return sensorID;
	}
	public ArrayList<Double> getSpeedData() {
		return speedData;
	}
	public void setSpeedData(ArrayList<Double> speedData) {
		this.speedData = speedData;
	}
	public ArrayList<Double> getDensityData() {
		return densityData;
	}
	public void setDensityData(ArrayList<Double> densityData) {
		this.densityData = densityData;
	}
	public ArrayList<Double> getFlowData() {
		return flowData;
	}
	public void setFlowData(ArrayList<Double> flowData) {
		this.flowData = flowData;
	}
	public void setSensorID(int sensorID) {
		this.sensorID = sensorID;
	}
	public CharSequence getSensorType() {
		return sensorType;
	}
	public void setSensorType(CharSequence charSequence) {
		this.sensorType = charSequence;
	}
	public String getSourceFileName() {
		return sourceFileName;
	}
	public void setSourceFileName(String sourceFileName) {
		this.sourceFileName = sourceFileName;
	}
	public String getSourceAddress() {
		return sourceAddress;
	}
	public void setSourceAddress(String sourceAddress) {
		this.sourceAddress = sourceAddress;
	}
	public double getHealthStatus() {
		return healthStatus;
	}
	public void setHealthStatus(double healthStatus) {
		this.healthStatus = healthStatus;
	}
	
	// methods
	public void addDatumToSpeed(double x){
		this.speedData.add(x);
	}
	
	public void addDatumToDensity(double x){
		this.densityData.add(x);
	}
	
	public void addDatumToFlow(double x){
		this.flowData.add(x);
	}
	
}