"""KIARA-RT AGN spectral table from qsosed (Kubota & Done 2018; pure-Python port by
Quera-Bofarull, PyPI 'qsosed', patched for astropy u.Ry and scipy simpson).

For each (log10 M_BH/Msun, log10 lambda = L/L_Edd, spin a*) the share of the bolometric
luminosity in each KIARA-RT photon group, with the spectrum cut at E_max (1 keV: harder
photons are left to the X-ray feedback), and the mean photon energy in each group.
lambda is qsosed's mdot (= L/L_Edd with qsosed's own efficiency). Points without a corona
solution (lowest lambda at high spin) take the values of the next lambda. The grid starts at
lambda = 10^-1.5: below that qsosed's corona solution breaks down (the shares can add up to
more than 1), and Obsidian's black holes are in their ADAF state, which emits no RT photons.
Used by BlackHoles:rt_agn_table_path (KIARA-RT with Obsidian).

The table also carries the properties of the optional AGN hard photon group of KIARA-RT
(the last edge of KIARART:photon_groups_Hz repeated, a build with one group more): the energy
and number weighted photoionization cross sections of HI, HeI and HeII (Verner et al. 1996,
as in src/rt/KIARA/rt_interaction_cross_sections.h) and the mean photon energy over
54.4 eV - E_max of a reference spectrum, the sum of the spectra (each normalised to the same
bolometric luminosity) of a grid of typical black holes (HARD_REF_*).

Usage: make_kiara_rt_agn_table.py out.hdf5 [nproc]   (about 4000 models, ~3 s each)
       make_kiara_rt_agn_table.py --add-hard-group table.hdf5 [nproc]   (60 models: adds
       the hard group attributes to an existing table)"""
import sys, itertools, warnings, numpy as np, h5py
from multiprocessing import Pool
warnings.filterwarnings("ignore")
from qsosed.sed import SED

EV = 1.602176634e-12
H_EV = 4.135667696e-15  # eV s
GROUPS_HZ = [1.451e15, 3.288e15, 5.945e15, 13.157e15]  # KIARA_4 group lower edges
E_MAX_EV = 1000.
edges = [nu * H_EV for nu in GROUPS_HZ] + [E_MAX_EV]
logM = np.arange(5.0, 10.0001, 0.25)
loglam = np.round(np.arange(-1.5, 0.3001, 0.1), 2)  # below 10^-1.5 qsosed's corona breaks down; Obsidian is in its ADAF state there
spin = np.array([-0.998, -0.9, -0.5, 0.0, 0.3, 0.5, 0.7, 0.9, 0.95, 0.998])

# Reference spectrum of the AGN hard group: black holes that dominate the AGN photons
HARD_REF_LOGM = [6.5, 7.0, 7.5, 8.0, 8.5]
HARD_REF_LOGLAM = [-1.5, -1.0, -0.5, 0.0]
HARD_REF_SPIN = [0.0, 0.5, 0.9]
# Verner et al. 1996 fits via Rosdahl et al. 2013, as in rt_interaction_cross_sections.h:
# E_ion, E_0 (eV), sigma_0 (cm^2), P, y_a, y_w, y_0, y_1
VERNER = {'HI': (13.6, 6.886e-13 / EV, 5.475e-14, 2.963, 32.88, 0., 0., 0.),
          'HeI': (24.59, 2.181e-11 / EV, 9.492e-16, 3.188, 1.469, 2.039, 0.4434, 2.136),
          'HeII': (54.42, 2.756e-12 / EV, 1.369e-14, 2.963, 32.88, 0., 0., 0.)}

def sigma(E, species):
    Ei, E0, s0, P, ya, yw, y0, y1 = VERNER[species]
    x = E / E0 - y0; y = np.sqrt(x * x + y1 * y1)
    return np.where(E >= Ei, s0 * ((x - 1.) ** 2 + yw ** 2) * y ** (0.5 * P - 5.5)
                    * (1. + np.sqrt(y / ya)) ** (-P), 0.)

def spectrum(args):
    lm, ll, a = args
    s = SED(M=10 ** lm, mdot=10 ** ll, spin=abs(a), spin_sign=(1 if a >= 0 else -1), reprocessing=True)
    E = s.ENERGY_RANGE_KEV * 1e3; L = s.total_flux(1.0) * s.ENERGY_BIN_WIDTHS_KEV * 1e3
    return E, L / L.sum()

def hard_group_properties(nproc):
    """Cross sections (cm^2; HI, HeI, HeII) and mean photon energy (eV) of the AGN hard
    group: the reference spectrum over the band of the table's last group to E_MAX_EV"""
    grid = list(itertools.product(HARD_REF_LOGM, HARD_REF_LOGLAM, HARD_REF_SPIN))
    with Pool(nproc) as p:
        res = p.map(spectrum, grid)
    E = res[0][0]; L = sum(r[1] for r in res)
    m = (E >= edges[-2]) & (E < edges[-1]); e = E[m]; l = L[m]; n = l / e
    cse = [(l * sigma(e, sp)).sum() / l.sum() for sp in ('HI', 'HeI', 'HeII')]
    csn = [(n * sigma(e, sp)).sum() / n.sum() for sp in ('HI', 'HeI', 'HeII')]
    return np.array(cse), np.array(csn), l.sum() / n.sum()

def write_hard_group(f, nproc):
    cse, csn, emean = hard_group_properties(nproc)
    f.attrs['hard_group_cse_cm2'] = cse; f.attrs['hard_group_csn_cm2'] = csn
    f.attrs['hard_group_mean_photon_energy_eV'] = emean
    f.attrs['hard_group_band_eV'] = [edges[-2], edges[-1]]
    f.attrs['hard_group_reference'] = (f'sum of spectra normalised to L_bol: log10 M_BH {HARD_REF_LOGM}, '
                                       f'log10 lambda {HARD_REF_LOGLAM}, spin {HARD_REF_SPIN}')
    print(f"hard group: mean photon energy {emean:.1f} eV, cse (HI, HeI, HeII) {cse}, csn {csn}")

def run(args):
    lm, ll, a = args
    ng = len(GROUPS_HZ)
    try:
        s = SED(M=10 ** lm, mdot=10 ** ll, spin=abs(a), spin_sign=(1 if a >= 0 else -1), reprocessing=True)
        E = s.ENERGY_RANGE_KEV * 1e3; dE = s.ENERGY_BIN_WIDTHS_KEV * 1e3
        L = s.total_flux(1.0) * dE          # per bin, arbitrary normalisation
        Ltot = L.sum()                      # bolometric (whole qsosed range)
        frac = np.zeros(ng); emean = np.zeros(ng)
        for g in range(ng):
            m = (E >= edges[g]) & (E < edges[g + 1])
            eb = L[m].sum(); nb = (L[m] / E[m]).sum()
            frac[g] = eb / Ltot; emean[g] = eb / nb if nb > 0 else 0.
        return frac, emean, s.efficiency, True
    except Exception:
        return np.full(ng, np.nan), np.full(ng, np.nan), np.nan, False

if __name__ == "__main__":
    if sys.argv[1] == '--add-hard-group':
        nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 8
        with h5py.File(sys.argv[2], 'r+') as f:
            if not np.allclose(f.attrs['photon_groups_Hz'], GROUPS_HZ):
                sys.exit("the table was made for other photon groups")
            write_hard_group(f, nproc)
        sys.exit(0)
    out = sys.argv[1]; nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    grid = list(itertools.product(logM, loglam, spin))
    with Pool(nproc) as p:
        res = p.map(run, grid, chunksize=2)
    shape = (len(logM), len(loglam), len(spin))
    frac = np.array([r[0] for r in res]).reshape(shape + (len(GROUPS_HZ),))
    emean = np.array([r[1] for r in res]).reshape(shape + (len(GROUPS_HZ),))
    eta = np.array([r[2] for r in res]).reshape(shape)
    ok = np.array([r[3] for r in res]).reshape(shape)
    # fill points without a corona solution from the next lambda
    for j in range(len(loglam) - 2, -1, -1):
        bad = ~ok[:, j, :]
        frac[:, j, :][bad] = frac[:, j + 1, :][bad]; emean[:, j, :][bad] = emean[:, j + 1, :][bad]
        eta[:, j, :][bad] = eta[:, j + 1, :][bad]; ok[:, j, :] |= ok[:, j + 1, :]
    with h5py.File(out, 'w') as f:
        f['log10_M_BH_Msun'] = logM; f['log10_eddington_ratio'] = loglam; f['spin'] = spin
        f['luminosity_fraction'] = frac; f['mean_photon_energy_eV'] = emean
        f['qsosed_efficiency'] = eta; f['valid'] = ok
        f.attrs['photon_groups_Hz'] = GROUPS_HZ; f.attrs['group_edges_eV'] = edges
        f.attrs['E_max_eV'] = E_MAX_EV
        f.attrs['source'] = 'qsosed (Kubota & Done 2018), reprocessing=True, default corona/warm-skin parameters'
        f.attrs['axes'] = 'luminosity_fraction[log10_M_BH_Msun, log10_eddington_ratio, spin, group]'
        write_hard_group(f, nproc)
    print("wrote", out, "shape", frac.shape, "models without corona (filled):", (~np.array([r[3] for r in res])).sum())
