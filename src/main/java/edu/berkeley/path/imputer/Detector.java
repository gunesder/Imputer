package edu.berkeley.path.imputer;

import java.io.Serializable;
import java.util.ArrayList;

public class Detector implements Serializable{
	
	// fields
	private FDParameters fdParams = new FDParameters();
	private int sensorID = 0;
	private CharSequence sensorType;
	private String sourceFileName;
	private String sourceAddress;
	private double healthStatus;
	private ArrayList<Double> speedData = new ArrayList<Double>();
	private ArrayList<Double> densityData = new ArrayList<Double>();
	private ArrayList<Double> flowData = new ArrayList<Double>();
	private int freewayNumber;
	private CharSequence freewayDirection;
	private Double absolutePM;
	private Double latitude;
	private Double longitude;
	private Double detectorLength;
	private CharSequence detectorName;
	private int numberOfLanes;
	private int linkAssoc;
	
	// getters and setters
	public FDParameters getFdParams() {
		return fdParams;
	}
	public void setFdParams(FDParameters fdParams) {
		this.fdParams = fdParams;
	}
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
	public Double getAbsolutePM() {
		return absolutePM;
	}
	public void setAbsolutePM(Double absolutePM) {
		this.absolutePM = absolutePM;
	}
	public Double getLatitude() {
		return latitude;
	}
	public void setLatitude(Double latitude) {
		this.latitude = latitude;
	}
	public Double getLongitude() {
		return longitude;
	}
	public void setLongitude(Double longitude) {
		this.longitude = longitude;
	}
	public Double getDetectorLength() {
		return detectorLength;
	}
	public void setDetectorLength(Double detectorLength) {
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
		public Double[] getSpeedDataArray() {
			Double[] d = speedData.toArray(new Double[speedData.size()]);
			Double[] dd = new Double[d.length];
			for(int i=0;i<d.length;i++){  
				 dd[i] = d[i].doubleValue();  
			}  
			return dd;
		}		
	public void setSpeedData(ArrayList<Double> speedData) {
		this.speedData = speedData;
	}
	public ArrayList<Double> getDensityData() {
		return densityData;
	}
		public Double[] getDensityDataArray() {
			Double[] d = densityData.toArray(new Double[densityData.size()]);
			Double[] dd = new Double[d.length];
			for(int i=0;i<d.length;i++){  
				 dd[i] = d[i].doubleValue();  
			}  
			return dd;
		}
	public void setDensityData(ArrayList<Double> densityData) {
		this.densityData = densityData;
	}
	public ArrayList<Double> getFlowData() {
		return flowData;
	}
		public Double[] getFlowDataArray() {
			Double[] d = flowData.toArray(new Double[flowData.size()]);
			Double[] dd = new Double[d.length];
			for(int i=0;i<d.length;i++){  
				 dd[i] = d[i].doubleValue();  
			}  
			return dd;
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
	public Double getHealthStatus() {
		return healthStatus;
	}
	public void setHealthStatus(Double healthStatus) {
		this.healthStatus = healthStatus;
	}
	
	// methods
	public void addDatumToSpeed(Double x){
		this.speedData.add(x);
	}
	
	public void addDatumToDensity(Double x){
		this.densityData.add(x);
	}
	
	public void addDatumToFlow(Double x){
		this.flowData.add(x);
	}
	
	public Double[] getSmoothedDataArray(String str) {
		Double[] d = new Double[speedData.size()];
		if (str.equals("speed")){
			d = speedData.toArray(new Double[speedData.size()]);
		} else if (str.equals("flow")){
			d = flowData.toArray(new Double[flowData.size()]);
		} else if (str.equals("density")){
			d = densityData.toArray(new Double[densityData.size()]);
		} else {
			throw new IllegalArgumentException("Illegal data retrieval");
		}
		Double[] dd = new Double[d.length];
		for(int i=0;i<d.length;i++){  
			 dd[i] = d[i].doubleValue();  
		}  
		dd = MyUtilities.smoothFWM(dd, 11, 2);
		return dd;
	}		
	
}