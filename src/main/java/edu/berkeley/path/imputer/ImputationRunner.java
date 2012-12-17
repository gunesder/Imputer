package edu.berkeley.path.imputer;

import java.io.FileNotFoundException;

import javax.xml.bind.JAXBException;

import core.DatabaseException;
import edu.berkeley.path.beats.simulator.SiriusException;

/**
 * Main method to call the imputer
 *
 */
public class ImputationRunner 
{
	// hard coded input parameters:
	public static final org.joda.time.DateTime startTime = new org.joda.time.DateTime(
			  // YYYY, MM, DD, HH, MM, TIME ZONE
		         2010,  1,  2,  0,  0, org.joda.time.DateTimeZone.forID("America/Los_Angeles")
		    );
	public static final org.joda.time.Duration totalTime = org.joda.time.Duration.standardHours(24);
	public static final String inputFileName = "C:\\Users\\gsr04\\Workspace\\imputer\\NetworkAConfig_NE.xml";
	public static final String outputFileName = "C:\\Users\\gsr04\\Workspace\\imputer\\NetworkAConfig_NE_out.xml";
	
    public static void main ( String[] args ) throws JAXBException, FileNotFoundException, SiriusException, DatabaseException
    {
    	Imputer imp = new Imputer(inputFileName,outputFileName,startTime,totalTime);
    	//Scenario scenario = imp.readAndUnmarshallXML();
    	imp.createDetectorListFromMainScenario();
    	imp.readDataIntoDetectorListFromDatabase();
    	imp.marshallIntoXML(imp.getMainScenario());
    }
}