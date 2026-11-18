#!/usr/bin/env python3
"""
Build the KIARA-RT BPASS photon tables from a BPASS SED grid.

KIARA-RT (src/rt/KIARA/rt_stellar_emission_model.h) wants, for each photon
group, the cumulative number of photons emitted per solar mass of stars from
age 0 to age t, on a grid of 13 metallicities x 22 ages (0, then 1 to 100 Myr
in steps of 0.1 dex). Each group is stored as a pandas-style 2D dataset
/Table_<name>/block0_values [13, 22].

The SED grid is the caesar/FSPS-format HDF5 file (e.g. BPASS_Chab100.hdf5):
  spectra [13*43, 100000]  L_nu in Lsun/Hz per Msun, metallicity-major
  wavelengths [100000]     1..100000 Angstrom
  ages [43]                log10(age/yr) = 6.0, 6.1, ..., 10.2
  logZ [13]

The method reproduces the original tables (built by Zhen Xiang) to < 0.5% for
HI, HeI and HeII:
  - photon rate Q(t) in a band is the trapezoidal integral of L_nu/(h nu) over
    the 1 Angstrom wavelength samples whose wavelength lies in [lmin, lmax];
  - N(0) = 0, N(1 Myr) = Q(1 Myr) * 1 Myr, and N is the trapezoidal time
    integral of Q(t) after that;
  - Lsun = 3.846e33 erg/s.
The original bands were HI 505-912 A, HeI 229-505 A and HeII <= 229 A.

By default the script copies the existing ionizing tables unchanged from
--ref (so the ionizing emission of runs is not altered) and adds a FUV
(Habing band, 6-13.6 eV, 912-2066 A) table. With no --ref, all tables are
computed from the SED grid. It also prints the photon-weighted mean energy of
each band and the fraction of the FUV photons in the Lyman-Werner band
(11.2-13.6 eV, 912-1107 A), integrated over 0-100 Myr.

Usage:
  make_kiara_rt_bpass_tables.py SED.hdf5 OUT.h5 [--ref OLD_TABLES.h5]
"""

import argparse

import h5py
import numpy as np

H_PLANCK = 6.62607015e-27  # erg s
C_LIGHT = 2.99792458e10  # cm/s
LSUN = 3.846e33  # erg/s, as in the original tables
MYR = 3.15576e13  # s
HC_EV_A = 12398.42  # h c in eV Angstrom
EV = 1.602176634e-12  # erg

N_Z = 13
N_AGE_SED = 43
METALLICITY_NAMES = [
    b"em5", b"em4", b"001", b"002", b"003", b"004", b"006", b"008",
    b"010", b"014", b"020", b"030", b"040",
]
# Table ages in Myr: 0, then 10^0.0 ... 10^2.0
T_TABLE = np.concatenate([[0.0], 10.0 ** np.round(np.arange(0.0, 2.01, 0.1), 1)])
AGE_NAMES = [("%.2f Myr" % t).encode() for t in T_TABLE]

# Bands as inclusive wavelength ranges [lmin, lmax] in Angstrom.
# The FUV band starts above 912 A, which belongs to the HI band.
BANDS = {
    "FUV": (913.0, 2066.0),
    "HI": (505.0, 912.0),
    "HeI": (229.0, 505.0),
    "HeII": (0.0, 229.0),
}
LW_BAND = (913.0, 1107.0)


def photon_rate(wl, spec, lmin, lmax):
    """Photons/s per Msun in [lmin, lmax] (Angstrom) for each SED row."""
    m = (wl >= lmin) & (wl <= lmax)
    nu = C_LIGHT / (wl[m] * 1e-8)
    integrand = spec[:, m] * LSUN / (H_PLANCK * nu)
    # wavelength increases -> nu decreases; integrate over increasing nu
    return np.trapezoid(integrand[:, ::-1], nu[::-1], axis=1)


def energy_rate(wl, spec, lmin, lmax):
    """erg/s per Msun in [lmin, lmax] (Angstrom) for each SED row."""
    m = (wl >= lmin) & (wl <= lmax)
    nu = C_LIGHT / (wl[m] * 1e-8)
    return np.trapezoid(spec[:, m][:, ::-1] * LSUN, nu[::-1], axis=1)


def cumulative(q):
    """Cumulative photons per Msun on T_TABLE from Q at log age 6.0..8.0."""
    t = 10.0 ** np.round(np.arange(6.0, 8.01, 0.1), 1) / 1e6  # Myr
    n = np.zeros(len(T_TABLE))
    n[1] = q[0] * 1.0 * MYR
    for k in range(len(t) - 1):
        n[k + 2] = n[k + 1] + 0.5 * (q[k] + q[k + 1]) * (t[k + 1] - t[k]) * MYR
    return n


def write_table(f, name, values):
    g = f.create_group("Table_" + name)
    g.attrs["pandas_type"] = np.bytes_("frame")
    g.attrs["pandas_version"] = np.bytes_("0.15.2")
    g.attrs["encoding"] = np.bytes_("UTF-8")
    g.create_dataset("axis0", data=np.array(AGE_NAMES))
    g.create_dataset("axis1", data=np.array(METALLICITY_NAMES))
    g.create_dataset("block0_items", data=np.array(AGE_NAMES))
    g.create_dataset("block0_values", data=values)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sed", help="BPASS SED grid (caesar/FSPS HDF5 format)")
    ap.add_argument("out", help="output table file")
    ap.add_argument("--ref", help="existing tables to copy the ionizing groups from")
    args = ap.parse_args()

    with h5py.File(args.sed, "r") as g:
        wl = g["wavelengths"][:]
        ages = g["ages"][:]
        if g["spectra"].shape[0] != N_Z * N_AGE_SED or len(ages) != N_AGE_SED:
            raise SystemExit("unexpected SED grid shape %s" % (g["spectra"].shape,))
        if not np.allclose(ages[:21], np.arange(6.0, 8.01, 0.1)):
            raise SystemExit("unexpected SED ages %s" % ages[:21])
        # only the first 21 ages (<= 100 Myr) are needed
        rows = [iz * N_AGE_SED + ia for iz in range(N_Z) for ia in range(21)]
        spec = g["spectra"][rows, :]

    tables = {}
    mean_energy = {}
    for name, (lmin, lmax) in BANDS.items():
        q = photon_rate(wl, spec, lmin, lmax).reshape(N_Z, 21)
        e = energy_rate(wl, spec, lmin, lmax).reshape(N_Z, 21)
        tables[name] = np.array([cumulative(q[iz]) for iz in range(N_Z)])
        etot = np.array([cumulative(e[iz]) for iz in range(N_Z)])
        mean_energy[name] = etot[:, -1] / np.maximum(tables[name][:, -1], 1e-300)
    q_lw = photon_rate(wl, spec, *LW_BAND).reshape(N_Z, 21)
    n_lw = np.array([cumulative(q_lw[iz]) for iz in range(N_Z)])
    f_lw = n_lw[:, -1] / tables["FUV"][:, -1]

    if args.ref:
        with h5py.File(args.ref, "r") as r:
            for name in ("HI", "HeI", "HeII"):
                old = r["Table_%s/block0_values" % name][:]
                ratio = tables[name][:, 1:] / old[:, 1:]
                print("%-4s rebuilt/ref: %.4f - %.4f (copying ref)"
                      % (name, ratio.min(), ratio.max()))
                tables[name] = old

    with h5py.File(args.out, "w") as f:
        for name in ("HI", "HeI", "HeII", "FUV"):
            write_table(f, name, tables[name])
        f.attrs["source_sed"] = np.bytes_(args.sed)
        f.attrs["FUV_band_Angstrom"] = BANDS["FUV"]
        f.attrs["LW_band_Angstrom"] = LW_BAND
        f.attrs["FUV_LW_photon_fraction_0_100Myr"] = f_lw
        for name in BANDS:
            f.attrs["mean_photon_energy_eV_0_100Myr_" + name] = mean_energy[name] / EV

    np.set_printoptions(precision=3, linewidth=160)
    print("Photon-weighted mean energy over 0-100 Myr [eV], by metallicity:")
    for name in BANDS:
        print("  %-4s" % name, mean_energy[name] / EV)
    print("LW fraction of FUV photons over 0-100 Myr:", f_lw)
    print("FUV photons per Msun to 100 Myr:", tables["FUV"][:, -1])
    print("wrote", args.out)


if __name__ == "__main__":
    main()
