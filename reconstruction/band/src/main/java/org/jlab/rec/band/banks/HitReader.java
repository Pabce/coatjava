package org.jlab.rec.band.banks;

import static org.jlab.rec.band.util.BandIdentifiers.pmtId;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import org.jlab.detector.banks.RawDataBank;

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.rec.band.constants.CalibrationConstantsLoader;
import org.jlab.rec.band.hit.BandHitCandidate;

/**
 * @author Efrain Segarra, Florian Hauenstein
 * Inspiration: Herbie
 */
public class HitReader {

	// this method retrieves the candidate hits ie the adc/tdc signal from all the pmt for a given event
	public static  ArrayList<BandHitCandidate> getBandCandidates(DataEvent event) {

		Map<Integer,Integer>	fadcInt 	= new HashMap<>();
		Map<Integer,Integer>	fadcAmpl    = new HashMap<>();
		Map<Integer,Float>	fadcTimes	= new HashMap<>();
		Map<Integer,Double>	tdcTimes	= new HashMap<>();
		Map<Integer,Integer> fadcIndex  = new HashMap<>();
		Map<Integer,Integer> tdcIndex  = new HashMap<>();

		if(event==null) return new ArrayList<>();

		// Grab trigger phase for TDC vs FADC matching
		double triggerPhase = getTriggerPhase(event);

		// Check if this is a MC file
		double TDCscale = getTDCScale(event);

		// Check that the file has the dgtz bank for BAND.	 
		if(event.hasBank("BAND::adc")==false || event.hasBank("BAND::tdc")==false) {
			return new ArrayList<>();
		} 

        // OLD:  directly accessing raw banks
		//DataBank bankADC = event.getBank("BAND::adc");
		//DataBank bankTDC = event.getBank("BAND::tdc");

        // NEW:  accessing raw banks via RawDataBank
        RawDataBank bankADC = new RawDataBank("BAND::adc");
        RawDataBank bankTDC = new RawDataBank("BAND::tdc");
        bankADC.read(event);
        bankTDC.read(event);

		ArrayList<BandHitCandidate> candidates = new ArrayList<>();

		int nadc 	= bankADC.rows();   // number of adc values
		int ntdc 	= bankTDC.rows();   // number of tdc values

		// Loop over the event and add all the FADC information to arrays for processing
		for( int i = 0 ; i < nadc ; i++ ){ 
			int s	= bankADC.getByte("sector",i);  // one of the 5 sectors
			int l	= bankADC.getByte("layer",i);  // one of the 6 layers
			int c 	= bankADC.getShort("component", i);
			int o	= bankADC.trueOrder(i);

			int adc = bankADC.getInt("ADC",i); 
			int ampl= bankADC.getInt("amplitude", i);
			float ftdc = bankADC.getFloat("time",i);
			int ped = bankADC.getShort("ped",i);
			ampl = ampl - ped;
	
			if( adc <= 0 || ftdc <= 0 ) continue;

			int key = pmtId(s, l, c, o);

			// Check if this PMT has been stored before, and if it has, then
			// replace it only if it has a larger ADC value. Otherwise, add to map
			//NOTE F.H. Feb 17 2020: Should we do some update to also check Amplitude here?
			if (fadcInt.containsKey(key) ) { 
				if( fadcInt.get(key) < adc ) {
					fadcInt.put(key, adc);
					fadcAmpl.put(key,  ampl);
					fadcTimes.put(key, ftdc);
					fadcIndex.put(key, bankADC.trueIndex(i));
				}
			}
			else {
				fadcInt.put(key, adc);
				fadcAmpl.put(key, ampl);
				fadcTimes.put(key, ftdc);
				fadcIndex.put(key, bankADC.trueIndex(i));
			}	
		}

		// Now loop over the TDC information and try to find best match
		for( int j = 0 ; j < ntdc ; j++ ){
			int s = bankTDC.getByte("sector", j);
			int l = bankTDC.getByte("layer", j);
			int c = bankTDC.getShort("component", j);
			int o = bankTDC.trueOrder(j);

			int key = pmtId(s, l, c, o-2);

			double tdc = ((double)bankTDC.getInt("TDC",j) * 0.02345 / TDCscale );

			if( tdc <= 0 ) continue;

			// Make sure that we have FADC information for this PMT -- if not, skip it
			if( !fadcInt.containsKey(key) || !fadcTimes.containsKey(key) ) continue;

			// If we have already stored a TDC for this PMT, need to compare with
			// FADC time and take the smallest tdiff one.
			if (tdcTimes.containsKey(key) ) { 
				double thisDiff = fadcTimes.get(key) - tdc;
				double prevDiff = fadcTimes.get(key) - tdcTimes.get(key);

				if( Math.abs(thisDiff) < Math.abs(prevDiff) ) {
					tdcTimes.put(key, tdc);
					tdcIndex.put(key, bankTDC.trueIndex(j));
				}
			}
			else {
				tdcTimes.put(key, tdc);
				tdcIndex.put(key, bankTDC.trueIndex(j));
			}
		}


		// Now we can loop over keys in our map and just ask that each map has that key in order
		// to save it as a hit candidate. The FADC is likely to have LESS info than the TDC, so 
		// we will loop over the keys in FADC to save time (since we are looking for match of all)
		for (Integer keys : fadcTimes.keySet()){
			if( !tdcTimes.containsKey(keys) || !fadcInt.containsKey(keys) ) continue;
			String id = Integer.toString(keys);
			int sector 		= Integer.parseInt( id.substring(0,1) );
			int layer 		= Integer.parseInt( id.substring(1,2) );
			int component 		= Integer.parseInt( id.substring(2,3) );
			int order		= Integer.parseInt( id.substring(3,4) );

			int adc = fadcInt.get(keys);
			int ampl= fadcAmpl.get(keys);
			float ftdc = fadcTimes.get(keys);
			double tdc = tdcTimes.get(keys);
			int indexadc = fadcIndex.get(keys);
			int indextdc = tdcIndex.get(keys);

			BandHitCandidate newHit = new BandHitCandidate( sector,layer,component,order,
					adc, ampl,tdc, ftdc ,triggerPhase, indexadc, indextdc);
			candidates.add(newHit);
		}
		return  candidates;  
	}	

    /** C++ reference policy: retain exactly one accepted ADC and one TDC per PMT.
     * Counts are taken after ADC cuts, before pairing. Raw row indices are retained.
     * Tagged overlay orders are not nominal C++ input and are rejected explicitly.
     */
    public static ArrayList<BandHitCandidate> getCppCandidates(DataEvent event) {
        ArrayList<BandHitCandidate> result = new ArrayList<>();
        if (!event.hasBank("BAND::adc") || !event.hasBank("BAND::tdc")) return result;
        boolean mc = event.hasBank("MC::Event");
        DataBank adc = event.getBank("BAND::adc"), tdc = event.getBank("BAND::tdc");
        Map<Integer,ArrayList<Integer>> aa = new java.util.TreeMap<>(), tt = new java.util.TreeMap<>();
        for (int i=0; i<adc.rows(); i++) {
            int s=adc.getByte("sector",i), l=adc.getByte("layer",i), c=adc.getShort("component",i), o=adc.getByte("order",i);
            if (!physicalChannel(s,l,c) || o<0 || o>1 || (mc && l==6 && o==1)) continue;
            int a=adc.getInt("ADC",i), rawAmp=adc.getInt("amplitude",i), amp=rawAmp-adc.getShort("ped",i);
            float time=adc.getFloat("time",i);
            if (a==0 || time==0 || !Float.isFinite(time)) continue;
            if (mc ? a/1E4<0.5 : rawAmp>=4095 || amp<250) continue;
            aa.computeIfAbsent(pmtId(s, l, c, o),k->new ArrayList<>()).add(i);
        }
        for (int i=0; i<tdc.rows(); i++) {
            int s=tdc.getByte("sector",i), l=tdc.getByte("layer",i), c=tdc.getShort("component",i), o=tdc.getByte("order",i)-2;
            if (!physicalChannel(s,l,c) || o<0 || o>1) continue;
            int key=pmtId(s, l, c, o);
            if (aa.containsKey(key)) tt.computeIfAbsent(key,k->new ArrayList<>()).add(i);
        }
        double phase=0;
        if (!mc && event.hasBank("RUN::config")) {
            long stamp=event.getBank("RUN::config").getLong("timestamp",0);
            if (stamp!=-1) phase=4.0*((stamp+3)%6);
        }
        for (Map.Entry<Integer,ArrayList<Integer>> entry:aa.entrySet()) {
            ArrayList<Integer> times=tt.get(entry.getKey());
            if (entry.getValue().size()!=1 || times==null || times.size()!=1) continue;
            int a=entry.getValue().get(0), t=times.get(0), raw=tdc.getInt("TDC",t);
            // C++ integer compound assignment truncates MC channels before ns conversion.
            if (mc) raw=(int)(raw/1E4);
            result.add(new BandHitCandidate(adc.getByte("sector",a),adc.getByte("layer",a),adc.getShort("component",a),
                adc.getByte("order",a),adc.getInt("ADC",a),adc.getInt("amplitude",a)-adc.getShort("ped",a),
                raw*0.02345,adc.getFloat("time",a),phase,a,t));
        }
        return result;
    }

    public static boolean physicalChannel(int sector, int layer, int component) {
        if (sector<1 || sector>5 || layer<1 || layer>6) return false;
        int[] count=layer==5 ? new int[]{3,7,5,5,0} : new int[]{3,7,6,6,2};
        return component>=1 && component<=count[sector-1];
    }

	private static double getTriggerPhase( DataEvent ev) {
		double tPh = 0.;
		if(ev.hasBank("RUN::config")) {
			DataBank  bank = ev.getBank("RUN::config");                       
			long timeStamp = bank.getLong("timestamp", 0);
			if( CalibrationConstantsLoader.JITTER_CYCLES > 0 && timeStamp >= 0 ) {
				tPh = CalibrationConstantsLoader.JITTER_PERIOD *
					( ( timeStamp + CalibrationConstantsLoader.JITTER_PHASE ) % 
					  CalibrationConstantsLoader.JITTER_CYCLES);
			}
		}
		return tPh;
	}

	private static double getTDCScale( DataEvent ev ){
		double tdcScale = 1;
		if(ev.hasBank("MC::Event")) tdcScale = 1E4;
		return tdcScale;
	}
}
