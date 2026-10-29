package org.jlab.rec.band.constants;

import static org.jlab.rec.band.util.BandIdentifiers.barId;

import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;
import org.jlab.rec.band.banks.HitReader;
import org.jlab.utils.groups.IndexedTable;

/** Immutable, run-scoped C++ calibration contract. Tables must be installed as a
 * coherent variation before enabling bandReconstruction=cpp in production.
 * Geometry z is relative to the ideal target; motherZ is applied exactly once.
 */
public final class BandCalibration {
    public static final String PREFIX = "/calibration/band/reconstruction_v2/";
    public static final java.util.List<String> TABLES = java.util.List.of(
            PREFIX+"config", PREFIX+"bars", PREFIX+"time_walk_left", PREFIX+"time_walk_right");
    public final boolean simulation;
    public final double motherZ, attenuation;
    public final int hitCut, laserCut;
    private final Map<Integer,Bar> bars;

    public static final class Bar {
        public final int id, status;
        public final double length, x, y, z, tdcOffset, fadcOffset, tdcVelocity, fadcVelocity;
        public final double tdcPaddle, fadcPaddle, tdcLayer, fadcLayer, tdcGlobal, fadcGlobal, adcPerMeV;
        private final double[] left, right;
        private Bar(IndexedTable b, IndexedTable l, IndexedTable r, int s, int layer, int c) {
            int[] idx={s,layer,c}; id=barId(s, layer, c);
            status=integer(b,"status",idx);
            length=value(b,"length",idx); x=value(b,"x",idx); y=value(b,"y",idx); z=value(b,"z",idx);
            tdcOffset=value(b,"tdc_offset",idx); fadcOffset=value(b,"fadc_offset",idx);
            tdcVelocity=value(b,"tdc_velocity",idx); fadcVelocity=value(b,"fadc_velocity",idx);
            tdcPaddle=value(b,"tdc_paddle",idx); fadcPaddle=value(b,"fadc_paddle",idx);
            tdcLayer=value(b,"tdc_layer",idx); fadcLayer=value(b,"fadc_layer",idx);
            tdcGlobal=value(b,"tdc_global",idx); fadcGlobal=value(b,"fadc_global",idx);
            adcPerMeV=value(b,"adc_per_mev",idx);
            left=new double[24];right=new double[24];
            for(int i=0;i<24;i++){left[i]=value(l,"p"+i,idx);right[i]=value(r,"p"+i,idx);}
            if (length<=0 || (status!=0 && (adcPerMeV<=0 || (layer!=6 && (tdcVelocity<=0 || fadcVelocity<=0)))))
                throw new IllegalArgumentException("Invalid active BAND calibration for bar "+id);
        }
        /** Four-term amplitude time walk; preserve C++ slots and operation order. */
        public double timeWalk(int side, double amplitude) {
            double[] p=side==1 ? left:right;
            double f0=p[0]+p[1]/Math.pow(amplitude,p[2]);
            double f1=p[6]/amplitude+p[7]/(Math.exp((amplitude-p[8])/p[9])+p[10]);
            double f2=p[12]+p[13]*amplitude;
            double f3=p[20]==0 ? 0 : p[18]*Math.sin((amplitude-p[19])/p[20])+p[21]/Math.pow(amplitude,p[22]);
            double result=f0+f1+f2+f3;
            if (!Double.isFinite(result)) throw new IllegalArgumentException("Nonfinite BAND time walk for "+id);
            return result;
        }
    }

    public BandCalibration(Function<String,IndexedTable> source) {
        IndexedTable cfg=source.apply(TABLES.get(0)), b=source.apply(TABLES.get(1));
        IndexedTable l=source.apply(TABLES.get(2)), r=source.apply(TABLES.get(3));
        int mode=integer(cfg,"simulation",0,0,0);
        if (mode!=0 && mode!=1) throw new IllegalArgumentException("BAND simulation must be 0 or 1");
        simulation=mode==1;motherZ=value(cfg,"mother_z",0,0,0);attenuation=value(cfg,"attenuation",0,0,0);
        hitCut=integer(cfg,"hit_cut",0,0,0);laserCut=integer(cfg,"laser_cut",0,0,0);
        if(attenuation<=0 || hitCut<1 || laserCut<hitCut) throw new IllegalArgumentException("Invalid BAND config");
        Map<Integer,Bar> loaded=new HashMap<>();
        for(int s=1;s<=5;s++) for(int layer=1;layer<=6;layer++) for(int c=1;physical(s,layer,c);c++) {
            Bar bar=new Bar(b,l,r,s,layer,c);loaded.put(bar.id,bar);
        }
        bars=Map.copyOf(loaded);
    }
    private static boolean physical(int s,int l,int c){return HitReader.physicalChannel(s,l,c);}
    public Bar bar(int id) {
        Bar b=bars.get(id);
        if(b==null) throw new IllegalArgumentException("Unknown BAND bar "+id);
        return b;
    }
    private static double value(IndexedTable t,String field,int... idx) {
        if(t==null || !t.hasEntry(idx) || !t.getEntryMap().containsKey(field))
            throw new IllegalArgumentException("Missing BAND calibration field "+field+" at "+java.util.Arrays.toString(idx));
        double v=t.getDoubleValue(field,idx);
        if(!Double.isFinite(v)) throw new IllegalArgumentException("Nonfinite BAND calibration "+field);
        return v;
    }
    private static int integer(IndexedTable t,String field,int... idx) {
        double v=value(t,field,idx);
        if(v!=(int)v) throw new IllegalArgumentException("Nonintegral BAND calibration "+field);
        return (int)v;
    }
}
