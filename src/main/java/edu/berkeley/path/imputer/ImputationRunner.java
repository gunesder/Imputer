package edu.berkeley.path.imputer;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;

import javax.xml.bind.JAXBException;

import jxl.read.biff.BiffException;
import jxl.write.WriteException;
import jxl.write.biff.RowsExceededException;

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
	public static final org.joda.time.Duration totalTime = org.joda.time.Duration.standardHours(2);
//	public static final String inputFileName = System.getProperty("user.dir") + "\\NetworkAConfig_NE.xml";
//	public static final String outputFileName = System.getProperty("user.dir") + "\\NetworkAConfig_NE_out.xml";
	public static final String inputFileName = System.getProperty("user.dir") + "\\Aimsun\\ImputationCode_2282104.xml";
	public static final String outputFileName = System.getProperty("user.dir") + "\\Aimsun\\ImputationCode_2282104_out.xml";
	
    public static void main ( String[] args ) throws JAXBException, BeatsException, IOException, SQLException, BiffException, RowsExceededException, WriteException, IndexOutOfBoundsException
    {
    	Imputer imp = new Imputer(inputFileName,outputFileName,startTime,totalTime);
    	//Scenario scenario = imp.readAndUnmarshallXML();
    	imp.createDetectorListFromMainScenario();
    	imp.createNodeStructureFromMainScenario();
    	imp.createLinkStructureFromMainScenario();
    	imp.createMainlineLinkStructureFromMainScenario(1);
//    	imp.readDataIntoDetectorListFromDatabase();
    	imp.readDataIntoDetectorListFromSpreadSheet(System.getProperty("user.dir") + "\\Aimsun\\PeMS_5m_SANDAG_Imputer2_2282104.xls");
    	// imp.exportMainlineDataToText();
    	// imp.exportDetectors();
    	imp.calibrateFundamentalDiagrams();
//    	imp.readFundamentalDiagramsFromXML();
//    	imp.readFundamentalDiagramsFromXML_AIMSUN();
    	imp.createCellStructure();
    	imp.runImputation();
    	imp.splitMegaCells();
    	imp.writeDemandsAndSplitRatiosToSpreadSheet(System.getProperty("user.dir") + "\\Aimsun\\PeMS_5m_SANDAG_Imputer2_2282104.xls");
//    	imp.marshallIntoXML(imp.getMainScenario());
    }
}