package org.jlab.rec.band.util;

/** Shared BAND detector keys, distinct from the sequential IDs in output banks. */
public final class BandIdentifiers {
    private BandIdentifiers() {}

    /** Preserve the decimal bar encoding used by BAND calibration and C++ maps. */
    public static int barId(int sector, int layer, int component) {
        return sector * 100 + layer * 10 + component;
    }

    /**
     * PMT key with zero-based side (0 left, 1 right), as in the raw ADC order.
     * Callers using BandHitCandidate.GetSide() must subtract one first.
     * Channel/order validation remains the caller's responsibility.
     */
    public static int pmtId(int sector, int layer, int component, int side) {
        return barId(sector, layer, component) * 10 + side;
    }
}
