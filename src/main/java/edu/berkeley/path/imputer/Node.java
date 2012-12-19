package edu.berkeley.path.imputer;

public class Node {
	// fields
	private int nodeID;
	private String nodeType;
	private int[] inLinks;
	private int[] outLinks;
	
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
	public int[] getInLinks() {
		return inLinks;
	}
	public void setInLinks(int[] inLinks) {
		this.inLinks = inLinks;
	}
	public int[] getOutLinks() {
		return outLinks;
	}
	public void setOutLinks(int[] outLinks) {
		this.outLinks = outLinks;
	}

}
