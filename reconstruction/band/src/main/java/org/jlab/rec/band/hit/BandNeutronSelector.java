package org.jlab.rec.band.hit;

import static org.jlab.rec.band.util.BandIdentifiers.barId;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import org.jlab.rec.band.banks.HitReader;

/** BANDReco::storeHits eligibility and goodNeutronEvent blocking, evaluated
 * before float bank serialization. No electron, beta, or analysis skim cuts.
 */
public final class BandNeutronSelector {
    private BandNeutronSelector() {}
    public static final double MAIN_THRESHOLD = 2.0;
    public static final double VETO_THRESHOLD = (2.0 / 7.2) * 2.0;

    /** Indices always refer to the original reconstructed bar list, not the
     * compact eligible-hit list. blocker is the first blocker in detector order.
     */
    public record HitDecision(int row, int barId, boolean eligible, boolean survivor,
                              int blocker, double energy, double tof, double x, double y) {}
    /** decision: 0 rejected, 1 single survivor, 2 clustered pair. */
    public record Result(int decision, int eligibleCount, int survivorCount,
                         int leadRow, int leadId, List<HitDecision> hits) {
        public Result { hits = List.copyOf(hits); }
        public boolean accepted() { return decision != 0; }
    }

    public static Result select(List<BandHit> bars, double startTime) {
        if (!Double.isFinite(startTime)) throw new IllegalArgumentException("Nonfinite BAND start time");
        int n = bars.size();
        boolean[] eligible = new boolean[n], survivor = new boolean[n];
        double[] tof = new double[n];
        int[] blocker = new int[n], barIds = new int[n];
        ArrayList<Integer> ordered = new ArrayList<>();
        HashSet<Integer> identities = new HashSet<>();
        for (int i = 0; i < n; i++) {
            BandHit h = bars.get(i);
            barIds[i] = barId(h.GetSector(), h.GetLayer(), h.GetComponent());
            if (!HitReader.physicalChannel(h.GetSector(), h.GetLayer(), h.GetComponent())
                    || !identities.add(barIds[i])) throw new IllegalArgumentException("Invalid/duplicate BAND bar");
            if (h.GetEnergy() == null || !Double.isFinite(h.GetEnergy())
                    || !Double.isFinite(h.GetMeanTime_TDC()) || !Double.isFinite(h.GetX())
                    || !Double.isFinite(h.GetY())) throw new IllegalArgumentException("Invalid calibrated BAND hit");
            tof[i] = h.GetMeanTime_TDC() - startTime;
            if (!Double.isFinite(tof[i])) throw new IllegalArgumentException("Nonfinite BAND TOF");
            eligible[i] = h.GetEnergy() >= (h.GetLayer() == 6 ? VETO_THRESHOLD : MAIN_THRESHOLD);
            blocker[i] = -1;
            if (eligible[i]) ordered.add(i);
        }
        // The C++ caller iterates an ordered map, including for equal-time ties.
        ordered.sort(Comparator.comparingInt(i -> barIds[i]));
        ArrayList<Integer> survivors = new ArrayList<>();
        for (int i : ordered) {
            BandHit h = bars.get(i);
            if (h.GetLayer() == 6) continue;
            int delta = h.GetLayer() == 4 && (h.GetSector() == 5
                    || ((h.GetSector() == 3 || h.GetSector() == 4) && h.GetComponent() == 6)) ? 2 : 1;
            for (int j : ordered) {
                BandHit other = bars.get(j);
                if (other.GetLayer() - h.GetLayer() != delta) continue;
                boolean veto = other.GetLayer() == 6;
                if (Math.abs(h.GetY() - other.GetY()) > 8
                        || (!veto && Math.abs(other.GetX() - h.GetX()) > 15)
                        || Math.abs(tof[j] - tof[i]) > (veto ? 15 : 3)) continue;
                blocker[i] = j;
                break; // All eligible hits can block, even if themselves blocked.
            }
            survivor[i] = blocker[i] == -1;
            if (survivor[i]) survivors.add(i);
        }
        int lead = -1, decision = 0;
        if (survivors.size() == 1) {
            lead = survivors.get(0); decision = 1;
        } else if (survivors.size() == 2) {
            int a = survivors.get(0), b = survivors.get(1);
            BandHit h = bars.get(a), other = bars.get(b);
            if (h.GetLayer() == other.GetLayer() && Math.abs(h.GetY() - other.GetY()) <= 8
                    && Math.abs(h.GetX() - other.GetX()) <= 15 && Math.abs(tof[a] - tof[b]) <= 3) {
                lead = tof[a] < tof[b] ? a : b; decision = 2;
            }
        }
        ArrayList<HitDecision> details = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            BandHit h = bars.get(i);
            details.add(new HitDecision(i, barIds[i], eligible[i], survivor[i], blocker[i],
                    h.GetEnergy(), tof[i], h.GetX(), h.GetY()));
        }
        return new Result(decision, ordered.size(), survivors.size(), lead,
                lead < 0 ? -1 : barIds[lead], details);
    }
}
