package edu.berkeley.path.imputer;

import java.io.*;
import java.util.*; 

import javax.xml.*;
import javax.xml.bind.*;
import javax.xml.validation.SchemaFactory;

import org.joda.time.Interval;

import core.DatabaseException;

import edu.berkeley.path.beats.jaxb.*;
import edu.berkeley.path.beats.simulator.JaxbObjectFactory;
import edu.berkeley.path.beats.simulator.SiriusException;
import edu.berkeley.path.beats.util.ScenarioUtil;

import edu.berkeley.path.model_elements.PeMSAggregate;
import edu.berkeley.path.model_elements.PeMSStation;
import edu.berkeley.path.model_elements.PeMSStationAggregate;
import edu.berkeley.path.scenario_database_access.DBParams;
import edu.berkeley.path.scenario_database_access.PeMSStationAggregateReader;
import edu.berkeley.path.scenario_database_access.PeMSStationReader;

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
	private static HashMap<Integer, Detector> detectors = new HashMap<Integer, Detector>();
	private Interval timeInterval;
	
	// getters and setters
	public static HashMap<Integer, Node> getNodes() {
		return nodes;
	}
	public static void setNodes(HashMap<Integer, Node> nodes) {
		Imputer.nodes = nodes;
	}
	public static HashMap<Integer, Link> getLinks() {
		return links;
	}
	public static void setLinks(HashMap<Integer, Link> links) {
		Imputer.links = links;
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

	// constructors
	public Imputer(String inFileName, String outFileName, org.joda.time.DateTime startTime, org.joda.time.Duration totalTime) throws FileNotFoundException, JAXBException, SiriusException {
		inputFileName = inFileName;
		outputFileName = outFileName;
		mainScenario = this.readAndUnmarshallXML();
		timeInterval = new Interval(startTime, totalTime);
	}
	
	// methods
	/**
	 * Loads the XML schema as a resource
	 * @returns the schema
	 * @throws SiriusException
	 */
	public static javax.xml.validation.Schema getSchema() throws SiriusException {
		SchemaFactory factory = SchemaFactory.newInstance(XMLConstants.W3C_XML_SCHEMA_NS_URI);
		try {
			return factory.newSchema(ObjectFactory.class.getClassLoader().getResource("sirius.xsd"));
		} catch (org.xml.sax.SAXException exc) {
			throw new SiriusException(exc);
		}
	}
	
	/**
	 * Takes input XML file and unmarshalls it
	 * @throws JAXBException, FileNotFoundException, SiriusException
	 * @returns Scenario
	 */
	
	public edu.berkeley.path.beats.jaxb.Scenario readAndUnmarshallXML() throws JAXBException, FileNotFoundException, SiriusException {
				
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
	public void marshallIntoXML(Scenario scenarioToWrite) throws JAXBException, FileNotFoundException, SiriusException {
		
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
		int i,j;
		for (i=0; i<this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().size(); i++){
			int[] inputs = null, outputs = null;
			Node n = new Node();
			n.setNodeID(Integer.parseInt(this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().get(i).getId()));
			n.setNodeType(this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().get(i).getType());
			for (j=0; j<this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().get(i).getInputs().getInput().size(); j++){
				inputs[j] = Integer.parseInt(this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().get(i).getInputs().getInput().get(j).getLinkId());
				outputs[j] = Integer.parseInt(this.mainScenario.getNetworkList().getNetwork().get(0).getNodeList().getNode().get(i).getOutputs().getOutput().get(j).getLinkId());
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
		int i;
		for (i=0; i<this.mainScenario.getNetworkList().getNetwork().get(0).getLinkList().getLink().size(); i++){
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
		}
	}
	
	/**
	 * Reads the SensorList from mainScenario and populates the detectors hashmap
	 */
	public void createDetectorListFromMainScenario() {
		int i;
		String sensorIDString;
		for (i=0; i<this.mainScenario.getSensorList().getSensor().size(); i++){
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
	 * @throws DatabaseException 
	 */
	public void readDataIntoDetectorListFromDatabase() throws DatabaseException {
		PeMSStationAggregateReader stationAggregateReader = new PeMSStationAggregateReader(new DBParams());
		ArrayList<Long> vdsIDs = new ArrayList<Long>();
		
		for (int key: detectors.keySet()){
			vdsIDs.add((long) key);
		}
		List<PeMSStationAggregate> stationsAggregate = stationAggregateReader.read(this.timeInterval,vdsIDs,PeMSAggregate.AggregationLevel.PEMS_5MIN);
		
		// Read 5 minute data into the hashmap
		int i;
		for (i=0; i<stationsAggregate.size(); i++){
			// find the detector corresponding to the current ID in the data vector and fill the fields accordingly
			Detector d = detectors.get(stationsAggregate.get(i).getVdsId().intValue());
			d.addDatumToSpeed(stationsAggregate.get(i).getTotal().getAvgSpeed());
			d.addDatumToFlow(stationsAggregate.get(i).getTotal().getFlow());
			d.addDatumToDensity(stationsAggregate.get(i).getTotal().getAvgOccupancy());
			if(i<detectors.size()){
				d.setHealthStatus(stationsAggregate.get(i).getTotal().getObserved());
			}
		}
		
		// Read absolute detector info into the hashmap
		PeMSStationReader stationReader = new PeMSStationReader(new DBParams());
		for (int key: detectors.keySet()){
			PeMSStation station = stationReader.read((long) key);
			Detector d = detectors.get(key);
			d.setAbsolutePM(station.getAbsPostmile());
			d.setDetectorLength(station.getDetectorLength());
			d.setDetectorName(station.getDetectorName());		
			d.setFreewayDirection(station.getDirection());
			d.setFreewayNumber(station.getFwyNum());
			d.setLatitude(station.getLatitude());
			d.setLongitude(station.getLongitude());
			d.setNumberOfLanes(station.getLaneCount());
		}
				
	}
	
	/**
	 * Translates the link structure into the cell structure depending on healthy detector locations
	 */
	public void createCellStructure() {
		
	}
	
	

}
