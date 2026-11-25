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

Usage: make_kiara_rt_agn_table.py out.hdf5 [nproc]   (about 4000 models, ~3 s each)"""
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
    print("wrote", out, "shape", frac.shape, "models without corona (filled):", (~np.array([r[3] for r in res])).sum())
