package edu.berkeley.path.imputer;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;

import javax.xml.bind.JAXBException;

import jxl.read.biff.BiffException;

import core.DatabaseException;
import edu.berkeley.path.beats.simulator.BeatsException;

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
	public static final org.joda.time.Duration totalTime = org.joda.time.Duration.standardHours(1);
//	public static final String inputFileName = System.getProperty("user.dir") + "\\NetworkAConfig_NE.xml";
//	public static final String outputFileName = System.getProperty("user.dir") + "\\NetworkAConfig_NE_out.xml";
	public static final String inputFileName = System.getProperty("user.dir") + "\\NetworkAConfig_BaseCase_13.04.12.xml";
	public static final String outputFileName = System.getProperty("user.dir") + "\\NetworkAConfig_BaseCase_13.04.12_out.xml";
	
    public static void main ( String[] args ) throws JAXBException, BeatsException, IOException, SQLException, BiffException
    {
    	Imputer imp = new Imputer(inputFileName,outputFileName,startTime,totalTime);
    	//Scenario scenario = imp.readAndUnmarshallXML();
    	imp.createDetectorListFromMainScenario();
    	imp.createNodeStructureFromMainScenario();
    	imp.createLinkStructureFromMainScenario();
    	imp.createMainlineLinkStructureFromMainScenario();
//    	imp.readDataIntoDetectorListFromDatabase();
    	imp.readDataIntoDetectorListFromSpreadSheet();
    	// imp.exportMainlineDataToText();
    	// imp.exportDetectors();
//    	imp.calibrateFundamentalDiagrams();
    	// imp.readFundamentalDiagramsFromXML();
    	imp.readFundamentalDiagramsFromXML_AIMSUN();
    	imp.createCellStructure();
    	imp.runImputation();
    	imp.splitMegaCells();
    	imp.marshallIntoXML(imp.getMainScenario());
    }
}