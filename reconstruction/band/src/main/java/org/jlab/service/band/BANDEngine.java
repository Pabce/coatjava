package org.jlab.service.band;

import java.util.ArrayList;
import java.util.Arrays;
import org.jlab.clas.reco.ReconstructionEngine;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.rec.band.constants.CalibrationConstantsLoader;
import org.jlab.rec.band.constants.Parameters;
import org.jlab.rec.band.banks.HitReader;
import org.jlab.rec.band.banks.RecoBankWriter;
import org.jlab.rec.band.hit.BandHit;
import org.jlab.rec.band.hit.BandHitCandidate;
import org.jlab.rec.band.hit.BandHitFinder;


/** @author Florian Hauenstein, Efrain Segarra
 * Service to return reconstructed BAND Hits - the output is in Hipo format
 */

public class BANDEngine extends ReconstructionEngine {

	public BANDEngine() {
		super("BAND", "hauensteinsegarra", "1.0");
	}

	int Run = -1;
    private boolean cppReconstruction;
    private org.jlab.rec.band.constants.BandCalibration cppCalibration;
    private int cppRun = -1;
    private org.jlab.detector.calib.utils.ConstantsManager cppConstants;
    private String cppVariation="default", cppTimestamp="";
    private void resetCppConstants() {
        cppConstants=new org.jlab.detector.calib.utils.ConstantsManager(cppVariation);
        cppConstants.setTimeStamp(cppTimestamp);
        cppConstants.init(org.jlab.rec.band.constants.BandCalibration.TABLES);
        cppCalibration=null;
    }

    @Override public org.jlab.detector.calib.utils.ConstantsManager getConstantsManager() {
        return cppReconstruction ? cppConstants : super.getConstantsManager();
    }
    @Override public void setVariation(String variation) {
        if(cppReconstruction) { cppVariation=variation; resetCppConstants(); }
        else super.setVariation(variation);
    }
    @Override public void setTimeStamp(String timestamp) {
        if(cppReconstruction) { cppTimestamp=timestamp; resetCppConstants(); }
        else super.setTimeStamp(timestamp);
    }

	@Override
		public boolean processDataEventUser(DataEvent event) {
			if(cppReconstruction) {
                if(!event.hasBank("RUN::config")) return false;
                int run=event.getBank("RUN::config").getInt("run",0);
                if(run<=0) return true;
                if(cppCalibration==null || cppRun!=run) {
                    cppCalibration=new org.jlab.rec.band.constants.BandCalibration(
                        table->getConstantsManager().getConstants(run,table));
                    cppRun=run;
                }
                if(event.hasBank("MC::Event")!=cppCalibration.simulation)
                    throw new IllegalArgumentException("BAND data/MC calibration mismatch for run "+run);
                ArrayList<BandHitCandidate> pmts=HitReader.getCppCandidates(event);
                ArrayList<BandHit> bars=org.jlab.rec.band.hit.BandReconstruction.find(pmts,cppCalibration);
                RecoBankWriter.appendBANDBanks(event,pmts,bars);
                return true;
            }
            //System.out.println("**** NEW EVENT ****");
			// update calibration constants based on run number if changed
			setRunConditionsParameters(event);

			ArrayList<BandHitCandidate> candidates = new ArrayList<BandHitCandidate>();   
			ArrayList<BandHit> hits = new ArrayList<BandHit>();

			//1) Search for valid PMT hits based on FADC/TDC for each PMT
			candidates = HitReader.getBandCandidates(event)	;	
			// exit if candidates list is empty, neither BAND::rawhits nor BAND::hits is filled in this case
			if(candidates.size()==0 )
				return true;

			//2) Find the BAND bar hits from the candidates
			BandHitFinder hitFinder = new BandHitFinder();
			hits = hitFinder.findGoodHits(candidates);

			//3) Write candidates and hits to the banks. 
			RecoBankWriter.appendBANDBanks(event,candidates,hits);


			return true;
		}

	@Override
    public boolean init() {
            String mode=getEngineConfigString("bandReconstruction","legacy");
            if(!mode.equals("legacy") && !mode.equals("cpp"))
                throw new IllegalArgumentException("Unknown BAND reconstruction mode: "+mode);
            cppReconstruction=mode.equals("cpp");
            cppCalibration=null;cppRun=-1;
            if(cppReconstruction) {
                resetCppConstants();
                registerOutputBank("BAND::hits","BAND::rawhits","BAND::laser");
                return true;
            }
			
	    String[]  bandTables = new String[]{
				"/calibration/band/time_jitter",
				"/calibration/band/lr_offsets",
				"/calibration/band/effective_velocity",
				"/calibration/band/paddle_offsets",
				"/calibration/band/layer_offsets",
				"/calibration/band/paddle_offsets_tdc",
				"/calibration/band/layer_offsets_tdc",
				"/calibration/band/attenuation_lengths",
				"/calibration/band/time_walk_amp_left",
				"/calibration/band/time_walk_amp_right",
				"/calibration/band/global_offsets",
				"/calibration/band/cuts",
				"/calibration/band/energy_conversion"
				//"/calibration/band/time_walk_corr_left",
				//"/calibration/band/time_walk_corr_right",
        };
    
		requireConstants(Arrays.asList(bandTables));
        this.registerOutputBank("BAND::hits","BAND::rawhits","BAND::laser");
		return true;
	}    
        
    @Override
    public void detectorChanged(int runNumber) {}

	public void setRunConditionsParameters(DataEvent event) {
		if(event.hasBank("RUN::config")==false) {
			System.err.println("RUN CONDITIONS NOT READ!");
		}
		else {
			int newRun = Run;        

			DataBank bank = event.getBank("RUN::config");
			newRun = bank.getInt("run", 0);  
			// Load the constants
			//-------------------
			if(Run!=newRun) {
				CalibrationConstantsLoader.Load(newRun,this.getConstantsManager()); 
				Run = newRun;
				Parameters.CreateGeometry(); // loading BAND params
			}
		}

	}

	public static void main (String arg[]) {

		BANDEngine en = new BANDEngine();
		en.init();


		String input = "bandtest.hipo";

		HipoDataSource  reader = new HipoDataSource();
		reader.open(input);
		String outputFile="test.hipo";
		HipoDataSync  writer = new HipoDataSync();
		writer.open(outputFile);
		
		int nevent = 0;

		while(reader.hasEvent() && nevent<2) {
			DataEvent event = (DataEvent) reader.getNextEvent();
			//System.out.println("***********  NEXT EVENT ************");
			//event.show();
			//if (event.hasBank("band::adc") && event.hasBank("band::tdc")){
			//	event.getBank("band::adc").show();
			//	event.getBank("band::tdc").show();
			//}
			en.processDataEventUser(event);
			writer.writeEvent(event);
			//event.getBank("band::hits").show();
			nevent++;

		}		
		writer.close();


	}
}
