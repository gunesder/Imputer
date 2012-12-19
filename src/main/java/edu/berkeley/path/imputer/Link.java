package edu.berkeley.path.imputer;

import edu.berkeley.path.beats.jaxb.Scenario;

/** Link class to represent the network as a succession of links
 * a simpler representation with respect to the network under the scenario
 */

public class Link {
	
	// fields
	private int linkID;
	private Node upNode;
	private Node downNode;
	private int[] upLinks;
	private int[] downLinks;
	private double length;
	private int lanesML;
	private int lanesHOV;
	private boolean hasDetector;
	private Detector detectorML;
	private Detector detectorHOV;
	
	// getters and setters
	public int getLinkID() {
		return linkID;
	}
	public void setLinkID(int linkID) {
		this.linkID = linkID;
	}
	public Node getUpNode() {
		return upNode;
	}
	public void setUpNode(Node upNode) {
		this.upNode = upNode;
	}
	public Node getDownNode() {
		return downNode;
	}
	public void setDownNode(Node downNode) {
		this.downNode = downNode;
	}
	public int[] getUpLinks() {
		return upLinks;
	}
	public void setUpLinks(int[] upLinks) {
		this.upLinks = upLinks;
	}
	public int[] getDownLinks() {
		return downLinks;
	}
	public void setDownLinks(int[] downLinks) {
		this.downLinks = downLinks;
	}
	public double getLength() {
		return length;
	}
	public void setLength(double length) {
		this.length = length;
	}
	public int getLanesML() {
		return lanesML;
	}
	public void setLanesML(int lanesML) {
		this.lanesML = lanesML;
	}
	public int getLanesHOV() {
		return lanesHOV;
	}
	public void setLanesHOV(int lanesHOV) {
		this.lanesHOV = lanesHOV;
	}
	public boolean isHasDetector() {
		return hasDetector;
	}
	public void setHasDetector(boolean hasDetector) {
		this.hasDetector = hasDetector;
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
	
	// methods

}