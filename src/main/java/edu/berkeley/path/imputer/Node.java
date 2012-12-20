package edu.berkeley.path.imputer;

import java.util.ArrayList;

public class Node {
	// fields
	private int nodeID;
	private String nodeType;
	private ArrayList<Integer> inLinks = new ArrayList<Integer>();
	private ArrayList<Integer> outLinks = new ArrayList<Integer>();
	
	// getters and setters
	public int getNodeID() {
		return nodeID;
	}
	public void setNodeID(int nodeID) {
		this.nodeID = nodeID;
	}
	public String getNodeType() {
		return nodeType;
	}
	public void setNodeType(String nodeType) {
		this.nodeType = nodeType;
	}
	public ArrayList<Integer> getInLinks() {
		return inLinks;
	}
	public void setInLinks(ArrayList<Integer> inLinks) {
		this.inLinks = inLinks;
	}
	public ArrayList<Integer> getOutLinks() {
		return outLinks;
	}
	public void setOutLinks(ArrayList<Integer> outLinks) {
		this.outLinks = outLinks;
	}

}
