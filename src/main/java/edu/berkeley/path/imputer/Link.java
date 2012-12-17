package edu.berkeley.path.imputer;

import edu.berkeley.path.beats.jaxb.Scenario;

/** Link class to represent the network as a succession of links
 * a simpler representation with respect to the network under the scenario
 */

public class Link {
	
	// fields
	private int linkID;
	private int upNode;
	private int downNode;
	private int[] upLinks;
	private int[] downLinks;
	private double length;
	private int lanesML;
	private int lanesHOV;
	private boolean hasDetector;
	private int detectorML;
	private int detectorHOV;
	
	// getters and setters
	public int getLinkID() {
		return linkID;
	}
	public void setLinkID(int linkID) {
		this.linkID = linkID;
	}
	public int getUpNode() {
		return upNode;
	}
	public void setUpNode(int upNode) {
		this.upNode = upNode;
	}
	public int getDownNode() {
		return downNode;
	}
	public void setDownNode(int downNode) {
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
	public int getDetectorML() {
		return detectorML;
	}
	public void setDetectorML(int detectorML) {
		this.detectorML = detectorML;
	}
	public int getDetectorHOV() {
		return detectorHOV;
	}
	public void setDetectorHOV(int detectorHOV) {
		this.detectorHOV = detectorHOV;
	}	
	
	// methods

}