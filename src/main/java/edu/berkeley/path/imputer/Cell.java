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
	
	
	
	
}
