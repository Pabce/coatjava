package org.jlab.rec.band.hit;

import java.util.ArrayList;
import java.util.Map;
import java.util.TreeMap;
import org.jlab.rec.band.constants.BandCalibration;
import org.jlab.rec.band.constants.Parameters;

/** Calibrated bars corresponding to BANDReco::createBars. Neutron eligibility
 * and blocking belong to a later stage, not to reconstruction or hit status.
 */
public final class BandReconstruction {
    private BandReconstruction() {}

    public static ArrayList<BandHit> find(ArrayList<BandHitCandidate> candidates, BandCalibration calibration) {
        Map<Integer,Integer> indices=new TreeMap<>();
        for(int i=0;i<candidates.size();i++) {
            BandHitCandidate p=candidates.get(i);
            int id=(p.GetSector()*100+p.GetLayer()*10+p.GetComponent())*10+p.GetSide()-1;
            if(indices.put(id,i)!=null) throw new IllegalArgumentException("Duplicate BAND PMT "+id);
        }
        Map<Integer,BandHit> bars=new TreeMap<>();
        for(Map.Entry<Integer,Integer> entry:indices.entrySet()) {
            int leftIndex=entry.getValue(), id=entry.getKey()/10;
            BandHitCandidate left=candidates.get(leftIndex);
            boolean veto=left.GetLayer()==6;
            if(!veto && left.GetSide()!=1) continue;
            Integer rightIndex=veto ? Integer.valueOf(leftIndex):indices.get(id*10+1);
            if(rightIndex==null) continue;
            BandHitCandidate right=candidates.get(rightIndex);
            BandCalibration.Bar b=calibration.bar(id);
            if(b.status==0) continue;
            double tl=left.GetTimeCorr()-(calibration.simulation ? 0:b.timeWalk(left.GetSide(),left.GetAmpl()));
            double tr=right.GetTimeCorr()-(calibration.simulation ? 0:b.timeWalk(right.GetSide(),right.GetAmpl()));
            if(!valid(left,tl) || !valid(right,tr)) continue;
            // Promote FADC floats before arithmetic, as in the C++ double fields.
            double fl=left.GetFtdc(), fr=right.GetFtdc();
            double dt=(tl-tr)-b.tdcOffset, df=(fl-fr)-b.fadcOffset;
            if(!veto && (Math.abs(dt)>b.length/b.tdcVelocity || Math.abs(df)>b.length/b.fadcVelocity)) continue;
            double time=0.5*(tl+tr-b.tdcOffset)-b.tdcPaddle-b.tdcLayer-b.tdcGlobal;
            double ftime=0.5*(fl+fr-b.fadcOffset)-b.fadcPaddle-b.fadcLayer-b.fadcGlobal;
            double x=veto ? b.x:b.x-(0.5*dt*b.tdcVelocity);
            // The reference uses global TDC x, the same sign on both sides, and
            // mu=1e10 cm. Preserve even this tiny factor for numerical parity.
            double factor=Math.exp((b.length/2.-x)/calibration.attenuation);
            double al=left.GetAdc()*factor, ar=right.GetAdc()*factor;
            double energy=Math.sqrt(al*ar)/b.adcPerMeV;
            if(!Double.isFinite(time) || !Double.isFinite(ftime) || !Double.isFinite(x) || !Double.isFinite(energy))
                throw new IllegalArgumentException("Nonfinite reconstructed BAND bar "+id);
            BandHit hit=new BandHit();
            hit.SetSector(left.GetSector());hit.SetLayer(left.GetLayer());hit.SetComponent(left.GetComponent());
            hit.SetMeanTime_TDC(time);hit.SetMeanTime_FADC(ftime);hit.SetDiffTime_TDC(dt);hit.SetDiffTime_FADC(df);
            hit.SetTLeft_TDC(tl);hit.SetTRight_TDC(tr);hit.SetTLeft_FADC(fl);hit.SetTRight_FADC(fr);
            hit.SetAdcLeft(al);hit.SetAdcRight(ar);hit.SetEnergy(energy);
            hit.SetX(x);hit.SetY(b.y);hit.SetZ(b.z+calibration.motherZ);
            hit.SetUx(veto ? b.length/Math.sqrt(12):0.5*b.tdcVelocity*0.3);
            hit.SetUy(Parameters.thickness/Math.sqrt(12));hit.SetUz((veto ? 2:Parameters.thickness)/Math.sqrt(12));
            // A single-ended veto links twice to its actual PMT, like the C++ bar.
            hit.SetIndexLpmt(leftIndex);hit.SetIndexRpmt(rightIndex);
            bars.put(id,hit); // Last PMT wins if both veto sides occur, as in C++.
        }
        ArrayList<BandHit> result=new ArrayList<>(bars.values());
        boolean veto=result.stream().anyMatch(h->h.GetLayer()==6);
        return new BandHitFinder().advancedHitFinder(result,veto,calibration.hitCut,calibration.laserCut);
    }
    private static boolean valid(BandHitCandidate p,double corrected) {
        return p.GetAmpl()!=0 && p.GetAdc()!=0 && p.GetFtdc()!=0 && p.GetTimeCorr()!=0 && corrected!=0;
    }
}
