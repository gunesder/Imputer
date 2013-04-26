package edu.berkeley.path.imputer;

import java.io.*;
import java.sql.SQLException;
import java.util.*; 

import javax.xml.*;
import javax.xml.bind.*;
import javax.xml.validation.SchemaFactory;

import jxl.NumberCell;
import jxl.Workbook;
import jxl.read.biff.BiffException;

import org.joda.time.Interval;

import core.DatabaseException;
import core.oraDatabase;

import edu.berkeley.path.beats.jaxb.*;
import edu.berkeley.path.beats.simulator.BeatsException;
import edu.berkeley.path.beats.simulator.JaxbObjectFactory;

import edu.berkeley.path.model_objects.measurements.PeMSAggregate;
import edu.berkeley.path.model_objects.measurements.VDS;
import edu.berkeley.path.model_objects.measurements.PeMSStationAggregate;
//import edu.berkeley.path.model_database_access.DBParams;
import edu.berkeley.path.model_database_access.TestConfiguration;
import edu.berkeley.path.model_database_access.measurements.PeMSStationAggregateReader;
import edu.berkeley.path.model_database_access.measurements.VDSReader;

/**
 * Top level imputer class
 * XML in - XML out for the time being
 */
public class Imputer {
	
	// fields
	private String inputFileName;
	private String outputFileName;
	private Scenario mainScenario;
	private static HashMap<Integer, Node> nodes = new HashMap<Integer, Node>();
	private static HashMap<Integer, Link> links = new HashMap<Integer, Link>();
	private static LinkedList<Link> mainlineLinks = new LinkedList<Link>();
	private static LinkedList<Cell> cells = new LinkedList<Cell>();
	private static HashMap<Integer, Detector> detectors = new HashMap<Integer, Detector>();
	private Interval timeInterval;
	public static double totalTimeInHours;
	
	// getters and setters
	public static HashMap<Integer, Node> getNodes() {
		return nodes;
	}
	public static double getTotalTimeInHours() {
		return totalTimeInHours;
	}
	public static void setTotalTimeInHours(double totalTimeInHours) {
		Imputer.totalTimeInHours = totalTimeInHours;
	}
	public static void setNodes(HashMap<Integer, Node> nodes) {
		Imputer.nodes = nodes;
	}
	public static LinkedList<Link> getmainlineLinks() {
		return mainlineLinks;
	}
	public static void setmainlineLinks(LinkedList<Link> links) {
		Imputer.mainlineLinks = links;
	}
	public Interval getTimeInterval() {
		return timeInterval;
	}
	public void setTimeInterval(Interval timeInterval) {
		this.timeInterval = timeInterval;
	}
	public HashMap<Integer, Detector> getDetectors() {
		return detectors;
	}
	public static void setDetectors(HashMap<Integer, Detector> detectors) {
		Imputer.detectors = detectors;
	}
	public Scenario getMainScenario() {
		return mainScenario;
	}
	public void setMainScenario(Scenario mainScenario) {
		this.mainScenario = mainScenario;
	}
	public String getInputFileName() {
		return inputFileName;
	}
	public void setInputFileName(String inputFileName) {
		this.inputFileName = inputFileName;
	}
	public String getOutputFileName() {
		return outputFileName;
	}
	public void setOutputFileName(String outputFileName) {
		this.outputFileName = outputFileName;
	}
	public static HashMap<Integer, Link> getLinks() {
		return links;
	}
	public static void setLinks(HashMap<Integer, Link> links) {
		Imputer.links = links;
	}
	public static LinkedList<Link> getMainlineLinks() {
		return mainlineLinks;
	}
	public static void setMainlineLinks(LinkedList<Link> mainlineLinks) {
		Imputer.mainlineLinks = mainlineLinks;
	}
	public static LinkedList<Cell> getCells() {
		return cells;
	}
	public static void setCells(LinkedList<Cell> cells) {
		Imputer.cells = cells;
	}

	// constructors
	public Imputer(String inFileName, String outFileName, org.joda.time.DateTime startTime, org.joda.time.Duration totalTime) throws FileNotFoundException, JAXBException, BeatsException {
		inputFileName = inFileName;
		outputFileName = outFileName;
		mainScenario = this.readAndUnmarshallXML();
		timeInterval = new Interval(startTime, totalTime);
		totalTimeInHours = (timeInterval.getEndMillis()-timeInterval.getStartMillis())/(1000*60*60);
	}
	
	// methods
	/**
	 * Loads the XML schema as a resource
	 * @returns the schema
	 * @throws SiriusException
	 */
	public static javax.xml.validation.Schema getSchema() throws BeatsException {
		SchemaFactory factory = SchemaFactory.newInstance(XMLConstants.W3C_XML_SCHEMA_NS_URI);
		try {
			return factory.newSchema(ObjectFactory.class.getClassLoader().getResource("sirius.xsd"));
		} catch (org.xml.sax.SAXException exc) {
			throw new BeatsException(exc);
		}
	}
	
	/**
	 * Takes input XML file and unmarshalls it
	 * @throws JAXBException, FileNotFoundException, SiriusException
	 * @returns Scenario
	 */
	
	public edu.berkeley.path.beats.jaxb.Scenario readAndUnmarshallXML() throws JAXBException, FileNotFoundException, BeatsException {
				
		JAXBContext jaxbContext = JAXBContext.newInstance("edu.berkeley.path.beats.jaxb");
		Unmarshaller jaxbUnmarshaller = jaxbContext.createUnmarshaller();
		javax.xml.validation.Schema schema = this.getSchema();
		jaxbUnmarshaller.setSchema(schema);
		edu.berkeley.path.beats.simulator.ObjectFactory.setObjectFactory(jaxbUnmarshaller, new JaxbObjectFactory());
		Scenario scenario = (Scenario) jaxbUnmarshaller.unmarshal(new FileInputStream(this.inputFileName));   
		return scenario; 
		
	}
	
	/**
	 * Marshalls a scenario object and writes into output XML file
	 * @throws JAXBException, SiriusException
	 */
	public void marshallIntoXML(Scenario scenarioToWrite) throws JAXBException, FileNotFoundException, BeatsException {
		
		JAXBContext jaxbContext = JAXBContext.newInstance("edu.berkeley.path.beats.jaxb");
		Marshaller jaxbMarshaller = jaxbContext.createMarshaller();
		jaxbMarshaller.setSchema(this.getSchema());
		jaxbMarshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);
		jaxbMarshaller.marshal(scenarioToWrite, new File(this.outputFileName));   
				
	}
	
	/**
	 * Reads the network geometry from mainScenario and populates the nodes hashmap
	 */
	public void createNodeStructureFromMainScenario() {
		for (int i=0; i<this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().size(); i++){
			ArrayList<Integer> inputs = new ArrayList<Integer>();
			ArrayList<Integer> outputs = new ArrayList<Integer>();
			Node n = new Node();
			n.setNodeID(Integer.parseInt(this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().get(i).getId()));
			n.setNodeType(this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().get(i).getType());
			for (int j=0; j<this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().get(i).getInputs().getInput().size(); j++){
				inputs.add(Integer.parseInt(this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().get(i).getInputs().getInput().get(j).getLinkId()));
			}
			for (int j=0; j<this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().get(i).getOutputs().getOutput().size(); j++){
				outputs.add(Integer.parseInt(this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().get(i).getOutputs().getOutput().get(j).getLinkId()));
			}
			n.setInLinks(inputs);
			n.setOutLinks(outputs);
			nodes.put(n.getNodeID(), n);
		}
	}
	
	/**
	 * Reads the network geometry from mainScenario and populates the links hashmap
	 */
	public void createLinkStructureFromMainScenario() {
		for (int i=0; i<this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().size(); i++){
			Link l = new Link(); boolean hasDetector = false; Detector detectorML = new Detector();
			l.setLinkType(this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().get(i).getType());
			l.setLinkID(Integer.parseInt(this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().get(i).getId()));
			l.setUpNode(nodes.get(Integer.parseInt(this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().get(i).getBegin().getNodeId())));
			l.setDownNode(nodes.get(Integer.parseInt(this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().get(i).getEnd().getNodeId())));
			l.setUpLinks(l.getUpNode().getInLinks());
			l.setDownLinks(l.getDownNode().getOutLinks());
			l.setLength(this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().get(i).getLength().doubleValue());
			l.setLanesML(this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().get(i).getLanes().intValue());
			for(int key: detectors.keySet()){
				if (detectors.get(key).getLinkAssoc() == l.getLinkID()){
					hasDetector = true;
					// Detector dummyDetector = detectors.get(key);
					// if (dummyDetector.getSensorType().toString().equals(arg0) //REVISIT
					detectorML = detectors.get(key);
				}
			}
			l.setHasDetector(hasDetector);
			l.setDetectorML(detectorML);
			
			links.put(l.getLinkID(),l);
			
		}
				
	}
	
	/**
	 * Reads the network geometry from mainScenario and populates the mainlineLinks list
	 */
	public void createMainlineLinkStructureFromMainScenario() {
		for (int i=0; i<this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().size(); i++){
			// collect only mainline links in the links list
			if (this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().get(i).getType().equals("freeway")){
				Link l = new Link(); boolean hasDetector = false; Detector detectorML = new Detector();
				l.setLinkID(Integer.parseInt(this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().get(i).getId()));
				l.setUpNode(nodes.get(Integer.parseInt(this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().get(i).getBegin().getNodeId())));
				l.setDownNode(nodes.get(Integer.parseInt(this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().get(i).getEnd().getNodeId())));
				l.setUpLinks(l.getUpNode().getInLinks());
				l.setDownLinks(l.getDownNode().getOutLinks());
				l.setLength(this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().get(i).getLength().doubleValue());
				l.setLanesML(this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().get(i).getLanes().intValue());
				for(int key: detectors.keySet()){
					if (detectors.get(key).getLinkAssoc() == l.getLinkID()){
						hasDetector = true;
						detectorML = detectors.get(key);
					}
				}
				l.setHasDetector(hasDetector);
				l.setDetectorML(detectorML);
				mainlineLinks.add(l);
			}
		}
		// sort mainline links
		mainlineLinks = this.recursiveLinkSort(mainlineLinks);
		
	}
	
	/**
	 * Sorts the static LinkedList mainlineLinks in the order they appear on the freeway
	 */
	private LinkedList<Link> recursiveLinkSort(LinkedList<Link> links2) {
		
		if (links2.size() == 1){
			
			return links2;
			
		} else {
			
			boolean swapMade = false;
			ListIterator<Link> itr1 = links2.listIterator();
			while (itr1.hasNext()){
				Link temp = itr1.next();
				int tempindex = itr1.nextIndex();
				// if this loop makes any switches, set the flag to true
				if (links2.getFirst().getUpNode().getNodeID() == temp.getDownNode().getNodeID()){
					swapMade = true;
					links2.addFirst(temp);
					links2.remove(tempindex);								
					return this.recursiveLinkSort(links2);
				}
			}
			
			if(!swapMade){
			// assign last n-1 links to links3
				LinkedList<Link> links3 = new LinkedList<Link>();
				Link temp = links2.getFirst();
				links2.removeFirst();				
				links3 = this.recursiveLinkSort(links2);
				links3.addFirst(temp);
				return links3;
			} else {
				return links2;
			}
						
		}
		
	}
	
	/**
	 * Reads the SensorList from mainScenario and populates the detectors hashmap
	 */
	public void createDetectorListFromMainScenario() {
		String sensorIDString;
		for (int i=0; i<this.mainScenario.getSensorList().getSensor().size(); i++){
			Detector d = new Detector();
			sensorIDString = this.mainScenario.getSensorList().getSensor().get(i).getParameters().getParameter().get(7).getValue();
			d.setSensorID(Integer.parseInt(sensorIDString));
			d.setSensorType(this.mainScenario.getSensorList().getSensor().get(i).getType());
			d.setSourceAddress("n/a");
			d.setLinkAssoc(Integer.parseInt(this.mainScenario.getSensorList().getSensor().get(i).getLinkReference().getId()));
			detectors.put(d.getSensorID(), d);
		}
	}
	
	/**
	 * Reads the data from database and writes into detectors hashmap
	 * @throws SQLException 
	 */
	public void readDataIntoDetectorListFromDatabase() throws SQLException {
//		TestConfiguration.dbSetup();
		
		PeMSStationAggregateReader stationAggregateReader = new PeMSStationAggregateReader(oraDatabase.doConnect());
		ArrayList<Long> vdsIDs = new ArrayList<Long>();
		
		for (int key: detectors.keySet()){
			vdsIDs.add((long) key);
		}
		List<PeMSStationAggregate> stationsAggregate = stationAggregateReader.read(this.timeInterval,vdsIDs,PeMSAggregate.AggregationLevel.PEMS_5MIN);
		
		// Read absolute detector info into the hashmap
		VDSReader stationReader = new VDSReader(oraDatabase.doConnect());
		for (int key: detectors.keySet()){
			VDS station = stationReader.read((long) key);
			Detector d = detectors.get(key);
			d.setAbsolutePM(station.getAbsolutePostmile());
			d.setDetectorLength(station.getDetectorLength());
			d.setDetectorName(station.getDetectorName());		
			d.setFreewayDirection(station.getDirection());
			d.setFreewayNumber(station.getFreewayNum());
			d.setLatitude(station.getPosition().getPoint().get(0).getLat());
			d.setLongitude(station.getPosition().getPoint().get(0).getLng());
			d.setNumberOfLanes(station.getLaneCount());
		}
		
		// Read 5 minute data into the hashmap
		for (int i=0; i<stationsAggregate.size(); i++){
			// find the detector corresponding to the current ID in the data vector and fill the fields accordingly
			Detector d = detectors.get((int) stationsAggregate.get(i).getVdsId());
			d.addDatumToSpeed(stationsAggregate.get(i).getTotal().getAvgSpeed());
			d.addDatumToFlow(stationsAggregate.get(i).getTotal().getFlow()*12/d.getNumberOfLanes()); // to get the hourly rate at 5 minute granularity, multiply by 12
			d.addDatumToDensity(stationsAggregate.get(i).getTotal().getFlow()*12/stationsAggregate.get(i).getTotal().getAvgSpeed()/d.getNumberOfLanes());
			if(i<detectors.size()){
				d.setHealthStatus(stationsAggregate.get(i).getTotal().getObserved());
			}
		}
								
	}
	
	/**
	 * Reads the detector data from spreadsheet and writes into detectors hashmap
	 * The files should be in the following format and placed in the root directory of the imputer project folder (for example, see detOutMainlines_431.csv)
	 * 1) 5 minute data granularity is assumed 
	 * 2) The data should be sorted by alphabetical order of detector IDs and the data column should be chronologically sorted for each detector
	 * @throws IOException 
	 * @throws BiffException 
	 */
	public void readDataIntoDetectorListFromSpreadSheet() throws BiffException, IOException {
		
		String filename = System.getProperty("user.dir") + "\\detOutMainlines_431.csv";
		Workbook workbook = Workbook.getWorkbook(new File(filename));
		
		int rowIndex = 1; // start the index at 1 and increase by number of data points after each iteration
		// Read absolute detector info and 5 minute data into the hashmap (some fields not important for fake detectors, left blank or 0 for the time being)
		for (int key: detectors.keySet()){
			Detector d = detectors.get(key);
			NumberCell nc = (NumberCell) workbook.getSheet(0).getCell(rowIndex, 4);
			d.setAbsolutePM(nc.getValue());
			d.setDetectorLength(0.0);
			d.setDetectorName(workbook.getSheet(0).getCell(rowIndex, 1).toString());		
			d.setFreewayDirection("");
			d.setFreewayNumber(0);
			d.setLatitude(0.0);
			d.setLongitude(0.0);
			NumberCell nc1 = (NumberCell) workbook.getSheet(0).getCell(rowIndex, 8);
			Double temp = nc1.getValue();
			d.setNumberOfLanes(temp.intValue());
			for (int k=rowIndex; k<=rowIndex+totalTimeInHours*60/5; k++){
				NumberCell ncSpeed = (NumberCell) workbook.getSheet(0).getCell(k, 6);
				NumberCell ncFlow = (NumberCell) workbook.getSheet(0).getCell(k, 5);
				d.addDatumToSpeed(ncSpeed.getValue());
				d.addDatumToFlow(ncFlow.getValue()*12/d.getNumberOfLanes()); // to get the hourly rate at 5 minute granularity, multiply by 12
				d.addDatumToDensity(ncFlow.getValue()*12/ncSpeed.getValue()/d.getNumberOfLanes());
			}
			d.setHealthStatus(1.0);
			rowIndex += totalTimeInHours*60/5;
		}
								
	}
	
	public void exportMainlineDataToText() throws IOException{
		
		for (int key: detectors.keySet()){
			
			Double[] flow = MyUtilities.scaleVector(detectors.get(key).getFlowDataArray(),(double) detectors.get(key).getNumberOfLanes());
			Double[] speed = detectors.get(key).getSpeedDataArray();
			Double[] density = MyUtilities.scaleVector(detectors.get(key).getDensityDataArray(),(double) detectors.get(key).getNumberOfLanes());
			
			PrintWriter outFlow = new PrintWriter(new FileWriter(key + "_flw.txt"));
			PrintWriter outSpeed = new PrintWriter(new FileWriter(key + "_spd.txt"));
			PrintWriter outDensity = new PrintWriter(new FileWriter(key + "_dty.txt"));
			
			for (int i=0;i<flow.length;i++){
				
				outFlow.println(flow[i]);
				outSpeed.println(speed[i]);
				outDensity.println(density[i]);
				
			}
			
			outFlow.close();
			outSpeed.close();
			outDensity.close();
			
		}
		
	}
	
	public void exportDetectors() throws IOException { // doesn't work, remove when done
		
		for (int key: detectors.keySet()){
			
			Detector d = detectors.get(key);
			FileOutputStream fileOut = new FileOutputStream("Detector " + key + ".ser");
		    ObjectOutputStream out = new ObjectOutputStream(fileOut);
		    out.writeObject(d);
		    out.close();
		    fileOut.close();
		          
		}
	}
	
	/**
	 * Calibrates fundamental diagram parameters for detectors
	 */
	public void calibrateFundemantalDiagrams() {
		for (int key: detectors.keySet()){
			FDCalibrator fdCalib = new FDCalibrator();
			detectors.put(key, fdCalib.calibrateParameters(detectors.get(key)));
		}
		
	}
	
	/**
	 * Translates the link structure into the cell structure depending on healthy detector locations
	 */
	public void createCellStructure() {
		
		int i = 0;
				
		while (i < mainlineLinks.size()-1){
			
			if (mainlineLinks.get(i).isHasDetector() & mainlineLinks.get(i).getDetectorML().getHealthStatus() == 100){
				Cell c = new Cell((int) totalTimeInHours*60/5);
				c.addLink(mainlineLinks.get(i));
				c.setDetectorML(mainlineLinks.get(i).getDetectorML());
				c.setDetectorHOV(mainlineLinks.get(i).getDetectorHOV());
				while (!mainlineLinks.get(i+1).isHasDetector() & i < mainlineLinks.size()-2 | (mainlineLinks.get(i+1).isHasDetector() & mainlineLinks.get(i+1).getDetectorML().getHealthStatus() != 100)){
					c.addLink(mainlineLinks.get(i+1));
					i++;
				}
				// Onramps and Offramps in the Cell
				for (Link l:c.getLinks()){
					
					c.addToOnrampPerLink(l.getUpNode().getInLinks().size()-1);
					c.addToOfframpPerLink(l.getDownNode().getOutLinks().size()-1);
					
					for (int linkID:l.getUpNode().getInLinks()){
						if (links.get(linkID).getLinkType().equals("onramp")){
							if (links.get(linkID).getDetectorML().getFlowData().isEmpty() | links.get(linkID).getDetectorML().getHealthStatus() != 100){
								c.addToImputeOR(true);
								c.appendZeroColumnToMeasuredOnrampFlow();
							} else {
								c.appendColumnToMeasuredOnrampFlow(links.get(linkID).getDetectorML().getFlowDataArray());
								c.addToImputeOR(false);
							}
						}
					}
					for (int linkID:l.getDownNode().getOutLinks()){
						if (links.get(linkID).getLinkType().equals("offramp")){
							if (links.get(linkID).getDetectorML().getFlowData().isEmpty() | links.get(linkID).getDetectorML().getHealthStatus() != 100){
								c.addToImputeFR(true);
								c.appendZeroColumnToMeasuredOfframpFlow();
							} else {
								c.appendColumnToMeasuredOfframpFlow(links.get(linkID).getDetectorML().getFlowDataArray());
								c.addToImputeFR(false);
							}
						}
						
					}
				}
				
				cells.add(c);
			}
						
			i++;
					
		}	
		i = 0;
	}
	
	/**
	 * Method to invoke the imputation algorithm
	 */
	public void runImputation(){
		
		ImputationCoreAlgorithm imp_algo = new ImputationCoreAlgorithm(cells,detectors);
		imp_algo.run();
		
	}
	

}
