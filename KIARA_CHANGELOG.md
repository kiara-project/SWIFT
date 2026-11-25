# KIARA changelog

This file records what each KIARA merge into `master` changed: bug fixes, new or changed runtime parameters, behaviour changes, and restart compatibility. It is meant as the record for `master`, so you don't have to read the individual development-branch commits. Commit hashes are given for the details.

Newest first. Normally every change goes through `romeel_dev`. Add a section here with each `romeel_dev` → `master` PR, and list work not yet merged under "Unreleased".

Covered: every `romeel_dev` merge (PRs #32–#55), plus `rt_fuv` (#53), which was merged directly. Not covered: the `zoom-kiara-merge` merges (#45, #50) and the direct upstream SWIFT merge 9636cc20c (2026-11-12). Pre-KIARA upstream SWIFT history is in the upstream repository.

| PR | date | merge | theme |
|---|---|---|---|
| [—](#unreleased) | — | — | Unreleased changes on `romeel_dev` |
| [#55](#pr-55) | 2026-11-23 | ddf619d5c | crackle/feedback/BH bug fixes, wind recoupling, HeII self-shielding, halo T_vir |
| [#54](#pr-54) | 2026-11-18 | 23224c129 | MAGMA2 gradient limiter, direct neighbour counts |
| [#53](#pr-53) | 2026-11-16 | 5abd30c94 | KIARA RT FUV band, G0 from RT and photon conservation |
| [#52](#pr-52) | 2026-11-12 | 7c533e9e1 | Firehose momentum-exchange fixes |
| [#51](#pr-51) | 2026-11-10 | a76aaad84 | Size-resolved dust, WN07 SF-law knobs, AGN-wind diagnostics |
| [#48](#pr-48) | 2026-11-04 | 6a85aa342 | KIARA_25 example, unified SF efficiency, optional photometry |
| [#47](#pr-47) | 2026-11-04 | 72aa91149 | Built-in crackle, H2/ISRF and dust physics, Sigma_SFR winds, BH dynamical friction |
| [#46](#pr-46) | 2026-10-24 | f44218be8 | Jet-recoupling runaway heating fix |
| [#44](#pr-44) | 2026-10-21 | 28d4bdf9a | Obsidian BH spin evolution; H2 across subgrid ISM exit |
| [#43](#pr-43) | 2026-10-18 | f5aebdc4d | RT sub-cycling deadlock and photon-conserving transport |
| [#42](#pr-42) | 2026-10-16 | c2215d710 | KIARA RT transport overhaul and variable reduced speed of light |
| [#41](#pr-41) | 2026-10-12 | e6e4e63bd | KIARA review fixes, formatting and documentation |
| [#40](#pr-40) | 2026-10-12 | dc942619a | Physics audit: cooling, firehose, enrichment, star formation and Obsidian fixes |
| [#38](#pr-38) | 2026-10-08 | af317b618 | Upstream sync with Obsidian accretion/jet rework, FIRE G0 fit and FOFGalaxyIDs |
| [#37](#pr-37) | 2026-07-30 | 6308b70f2 | Startup smoothing-length depth fix |
| [#36](#pr-36) | 2026-07-30 | 2ae14470c | Obsidian AGN wind rework, SF-gas FoF galaxies, new FIRE G0 fit |
| [#35](#pr-35) | 2026-05-13 | c88b100dd | Cooling-shutoff follow-ups and BH repositioning by galaxy mass |
| [#34](#pr-34) | 2026-05-05 | 11056dfb3 | AGN cooling shutoff for all modes |
| [#33](#pr-33) | 2026-05-04 | 5e09e5436 | Upstream SWIFT sync |
| [#32](#pr-32) | 2026-04-26 | 1b3186832 | Initial KIARA model import: MAGMA2, KIARA cooling/chemistry/feedback, Obsidian BHs, KIARA-RT |

---

<a id="unreleased"></a>

## Unreleased (on `romeel_dev`)

- **Added `KIARA_CHANGELOG.md`** (667d0bd9e).
- **Resolved leftover merge-conflict markers in `configure.ac`** (b897bf792). They sat in the `--with-subgrid` help string, from PR #38. The help text now lists the options the `case` statement actually accepts.
- **Removed `ObsidianAGN:reposition_central_only`** (b897bf792). It has had no effect since PR #34, and a parameter file that still sets it now gets it reported as unused.
- **New option `ObsidianAGN:jet_velocity_cold_km_s`** (9efb4cf66), default 0 (off). The resolved jet speed is interpolated in log from `jet_velocity_km_s` at host gas fraction 0 to this value at gas fraction 1, with f_gas = M_gas/(M_gas+M*) from the FOF galaxy. Jet power is unchanged; gas-rich hosts get heavier, slower jets (more subgrid entrainment). Restart compatible.
- **New option `KIARAFeedback:recouple_escape_length`** (3cc710659), default 0 (launch kernel, as before). With 1, the escape-time recoupling (`recouple_escape_kernels`) measures crossing times over the initial firehose stream radius, which scales with the galaxy size, for both stellar and BH winds. Recalibrate N when switching. Restart compatible.
- **New option: cosmic rays in the subgrid ISM** (a89ed967f), `KIARACooling:cosmic_ray_ionization_rate_MW`, default 0 (off). This is the ionization rate per H atom (s⁻¹) at the Milky Way sSFR (`cosmic_ray_ssfr_MW_per_Gyr`, default 0.027); it scales linearly with the host galaxy sSFR, with an optional cap `cosmic_ray_max_factor` (in MW units). In crackle it ionizes HI and HeI (1.1× per He atom) and heats the neutral gas by `cosmic_ray_heating_eV` (default 10) per ionization. Before this, nothing heated the cold phase of metal-free gas, so the self-enrichment metals cooled it to T_CMB, x_e fell to ~1e-7 and H⁻ (the only primordial H2 route) shut off. Not restart compatible with earlier builds (new fields in the dumped `cooling_function_data`).
- **New option: PopIII stars** (fd00ff8a6), `KIARAFeedback:popIII_max_metallicity`, default −1 (off). This exposes chem5's existing PopIII channel (`zmax3`, previously hardcoded off) together with its IMF: `popIII_imf_slope` 1.35 (Salpeter), `popIII_mass_min_Msun` 20, `popIII_mass_max_Msun` 120, with the Z = 0 SN/HN yields. Stars above 50 M☉ collapse without ejecta, and there are no pair-instability SNe. Below the same threshold (gas + dust metallicity), the WadaNorman and lognormal SF laws use an H2 fraction of at least `KIARAStarFormation:popIII_H2_fraction` (default 1), so metal-free gas can form stars without self-enrichment; the stars it forms get the PopIII IMF. Not restart compatible with earlier builds (new fields in the dumped `feedback_props` and `star_formation`). (fd00ff8a6 aborted at startup any run whose parameter file did not set `popIII_max_metallicity`, because the SF and chem5 code read it with different defaults; fixed in e35585d09.)
- **New option `KIARAStarFormation:lowZ_KMT_metallicity`** (e37ea1e27), solar units, default 0 (off). With the WadaNorman and lognormal SF laws, gas below this metallicity uses at least the Krumholz & Gnedin (2011) H2 fraction evaluated at it (surface density ρ²/|∇ρ|, `clumping_factor_scaling`), as Simba's KMT does with Z′ floored at 0.01. Below ~1% solar, stars form in cold atomic gas (Krumholz 2012), so the crackle H2 fraction no longer gates their formation. Metal-free gas then forms stars above Σ ≈ 70 M☉/pc² (Z′ = 0.01, clumping 30) instead of needing a high self-enrichment floor. The KMT calculation is now a shared helper; SchmidtLaw + `H2_model: KMT` is unchanged. Not restart compatible with earlier builds (new field in the dumped `star_formation`).
- **New option `KIARACooling:HeII_reionization_redshift`** (5ff4deac5), default −1 (off), with `HeII_reionization_width` 0.5. The UVB HeII photoionization and photoheating rates are multiplied by ½[1 − tanh((z − z_HeII)/width)]. An optically thin uniform UVB (FG2011) otherwise ionizes and heats the IGM HeII by z ≈ 4.5 under `self_shielding_method: 4`, well before the observed HeII reionization (z ≈ 3–3.5). Not restart compatible with earlier builds (new fields in the dumped `cooling_function_data` and chemistry data).


### Upstream SWIFT
- Merged upstream master in 9cb92007f. A tighter upper bound on `dx_max_sort` stays valid when cells are re-sorted in only some directions (7c94f9ab1, with the new unit test `testSortShiftBound`). The multipole check of a single-particle cell skips higher-order terms (f3049a168). Restart compatible: the new fields are in `struct cell`, which is rebuilt on restart.

---

<a id="pr-55"></a>

## PR #55 — crackle/feedback/BH bug fixes, wind recoupling, HeII self-shielding, halo T_vir (2026-11-23, merge ddf619d5c)

These are bug fixes from the Nov 2026 KIARA calibration (s12.5n128 / s25n256) in crackle, KIARA feedback/firehose, Obsidian BHs and MAGMA2. The merge also adds options for wind recoupling, HeII self-shielding and the halo T_vir used by the BH feedback.

> ⚠️ Restart files from earlier builds are not compatible. New spart fields come from 796164d66 and new part fields from 08a8b7bdc.

### Bug fixes
- **Photoelectric heating was ~12× too weak** (af2f9179b). KIARA set `photoelectric_heating = 0`, so crackle used a rate that already includes ε·G0 inside the Wolfire+95 form, which applies them again. KIARA now uses mode 2 (γ_h = 1e-24 erg/s).
- **Negative hydro du/dt was dropped by crackle** (e8a725f9b). Adiabatic/PdV cooling of non-ISM gas was lost. The term is now added with no sign test, as Grackle does.
- **The H2 self-shielding column was 2× too large** (cf984b2ab). It counted H nuclei in H2, not H2 molecules (inherited from Grackle). The Sternberg+14 cap now subtracts an atomic layer on both sides.
- **SN dust destruction compounded on every crackle iteration** (d9f17e2f9). It is now applied once per step.
- **Non-ISM dust was sputtered twice** (808375582), once in crackle and again in `cooling_sputter_dust()`.
- **Decoupled and cooling-shutoff gas kept frozen ISM ionization states** (0f60b8d0b).
  - Heated to 1e5–1e6 K while still neutral, this gas held 37% of the box HI in s25n256 at z=0. On recoupling it radiated ~90% of u within ~2 Myr.
  - Its species now evolve at fixed internal energy, via a new crackle `fixed_energy` mode.
- **Existing winds could be re-kicked by nearby stars** (3c97387e4). A re-kick reset their decoupling timer and took launch mass away from star-forming gas.
- **The SFR-weighted launch weight depended on the unit system** (796164d66). It was m(W + SFR). Winds now launch from star-forming gas with weight mW·SFR, as in Simba, or from all eligible gas with mW if no neighbour forms stars.
- **The cold-wind temperature was converted with μ = 0.59** (92493e1f4). It now uses the particle's own μ, so winds no longer start at ≥2× the intended T.
- **The firehose energy caps were not energy-conserving** (6ac579a8d). Capped energy was lost after the partner had already exchanged its share. The caps are removed, since the exchange is already bounded.
- **The MAGMA2 velocity Hessian was contracted with the wrong index** (b580df959). The velocity gradient tensor was also stored gradient-first.
  - The Hessian fix changes the second-order reconstruction term.
  - KH growth (N=64) changes from 3.03 to 2.86. Gresho, Noh and Sod are unchanged.
- **BH repositioning could target past star-forming gas** (560c1d9f8). KIARA marks such gas with SFR < 0, and it may now be hot or diffuse. Repositioning now requires SFR > 0, with a fallback for quenched galaxies.
- **The BH merger radius was in units of h_BH/2.8 (~0.5 h_BH)** instead of the softening (1e14e4c91).

### New features / options
| parameter | default | what it does |
|---|---|---|
| `KIARAFeedback:recouple_escape_kernels` | 0 (off) | Recouples a non-jet wind still (or again) in ISM-density ambient gas after N launch kernel-crossing times (08a8b7bdc). |
| `KIARACooling:self_shielding_method` | 3 | New value 4: HI/HeI as in method 3, plus Rahmati-style HeII self-shielding, so the IGM reionizes HeII (e536c635d). |
| `ObsidianAGN:virial_temperature_model` | 0 | 1 = halo T_vir from host M* via inverse Moster+13, instead of the M_BH fit (0c0528e8c). Output as `VirialTemperatures`. |
| `ObsidianAGN:merger_radius_from_softening` | 1 | Merger radius in softening units; 0 restores h_BH (1e14e4c91). |
| `ObsidianAGN:reposition_min_sf_ngbs` | 1 | Below this many star-forming neighbours, any gas is a repositioning target (560c1d9f8). |
| `KIARAChemistry:firehose_log_energy` | 0 | Diagnostic log of large firehose energy exchanges (c7aa4e10c). |

- **Wind recoupling** (08a8b7bdc). The low-density recoupling test now uses the ambient density of coupled neighbours, not the SPH density, which stuck winds inflated for each other. With N=4 in s12.5n128, decoupled gas at ISM density drops ~10×.
- **RT** (e536c635d). With method 3, HeII no longer absorbs RT photons in the subgrid path. Before, it destroyed them without ionizing or heating anything.
- **SPHENIX with KIARA** (6f6ef7205). The KIARA hooks were ported from MAGMA2, so SPHENIX builds and runs with KIARA.

### Behaviour changes
- **BH mergers use the softening as their radius by default** (`merger_radius_from_softening: 1`).
- **BH repositioning targets only currently star-forming gas.**
- **FOF skips foreign cells whose particles were not received**, with a warning instead of an abort (a1fadfe53). This is a temporary workaround: s25n256 hits it at z≈0.1–1.1.
- **The firehose energy caps are gone**, so results change wherever they used to act.

### Other
- Documented the `adaf_jet_free_energy` assumption; comments only (5d93fa0b9).
- Recommended template settings: `self_shielding_method: 4`, `self_enrichment_metallicity: 1.e-4`, `recouple_escape_kernels: 4.`, `virial_temperature_model: 1`. `thermostat_factor`, `lowz_energy_ramp_*` and `sigma_crit_Msun_pc2` are not read by this code.

---

<a id="pr-54"></a>

## PR #54 — MAGMA2 gradient limiter, direct neighbour counts (2026-11-18, merge 23224c129)

One commit (88184ed4b) that adds an optional MAGMA2 limiter on oversized matrix gradients, lets stars and BHs set their neighbour number directly, and adds two calibration knobs (metal-line cooling factor, separate BH-suppression mass loading).

> ⚠️ Restart files from earlier builds are not compatible. (New `part` fields `force_interaction_count` and `gradient_limited_count`.)

### New features / options
| parameter | default | what it does |
|---|---|---|
| `SPH:max_gradient_ratio` | 0 (no limit) | MAGMA2: a pair falls back to the SPH kernel gradient when \|G_ij\| exceeds this many times the pair's mean SPH kernel gradient. Must be >= 0 |
| `Stars:target_neighbours` | 0 (unset) | Sets the stellar kernel from a neighbour number instead of `Stars:resolution_eta`. Set at most one; if neither is set, the hydro value is used |
| `BlackHoles:target_neighbours` | 0 (unset) | The same for BHs, as an alternative to `BlackHoles:resolution_eta` |
| `KIARACooling:metal_line_cooling_factor` | 1 | Multiplies the Cloudy metal cooling/heating term in crackle. 0 removes metal-line cooling while keeping the tables loaded |
| `ObsidianAGN:suppression_eta_normalization` | `KIARAFeedback:FIRE_eta_normalization` | Mass loading used by BH growth suppression (`suppress_growth` 4 or 5), so the winds can be varied at fixed BH suppression |

- New MAGMA2 snapshot fields `ForceInteractionCounts` and `GradientLimitedCounts`: cumulative counts of coupled force interactions and of those that hit the gradient limiter.

### Other
- `kiara_25.yml`: drops `dust_growth_tauref` and sets `dust_growth_densref` to 2.0e-24 (was 5.0e-24).

---

<a id="pr-53"></a>

## PR #53 — KIARA RT FUV band, G0 from RT and photon conservation (2026-11-16, merge 5abd30c94)

Adds a non-ionizing FUV (Habing-band) photon group to KIARA RT, uses it to set G0 and LW H2 dissociation in crackle, and makes the RT transport, stellar emission and subgrid-ISM absorption conserve photons. This came from the rt_fuv branch, merged straight into master rather than through romeel_dev, and also brought in a batch of upstream SWIFT commits.

> ⚠️ Restart files from earlier builds are not compatible.

### Bug fixes
- **Transport was not conservative** (15d41b9a3). An inactive neighbour got the interface area of every step of a smaller-step active particle; each interface is now counted once per own step. Flux was converted to energy with the scale factor at the time it was applied rather than at the exchange, so inactive receivers got Δa/a less than the sender lost (about 3% of injected energy by z = 9).
- **Radiation was lost or stuck in ineligible particles** (08e13c44b, f76bae273). Decoupled winds and cooling-off particles now hand their radiation to eligible neighbours over their next RT step, conservatively. Stars inject into all neighbours, and a star with no neighbours keeps its emission for the next step. A gas particle turning into a star passes its radiation (including pending transport) to the star. The fraction handed over per RT step is capped at 1, which stops the negative-energy clamp creating energy (a runaway in the FUV runs).
- **RT gradient loop dropped neighbours of large-h particles** (d5e4729be). The RT loop cut-offs used `hydro.h_max_active`, which misses RT-active particles on RT sub-cycles (e.g. decoupled winds). The RT loops now use `hydro.h_max`; the hydro loops are unchanged.
- **Stellar emission injected all at once** (4b8cfc91c). A star's whole-step emission went into its neighbours in one go and left the dense gas within ~kyr, so star-forming gas had G0 10^3–10^4 below the stellar field. Neighbours now hold it as `emission_pending` and add rate·dt of it each RT step. In a 6.25 Mpc/h test to z = 9, G0 summed over star-forming gas went from 2e-4 to 0.3–0.6 of the unattenuated stellar field, with 16% more run time.
- **Untransported emission in the final snapshot** (c0259d37d). The final step has no RT sub-cycles, so its emission sat next to young stars (FUV G0 1e3–1e4 too high). The star time step is now set to zero on the final step.
- **RT crackle instance ignored KIARACooling** (455c00997). ISRF H2 photodissociation and the Sternberg limit were off, `H2_self_shielding` was 4 instead of 3, and dust parameters were hard-coded. The new `cooling_set_crackle_dust_parameters()` is applied to the RT instance after `cooling_init()` and on restart. This changes the results of existing KIARA-RT runs.
- **Huge interface areas at the smoothing-length floor** (c11f72e4a). Areas reached ~1000× those of a sphere, so particles exported most of their radiation every step. Areas are now capped in the exchange too, conservatively, using min(s_i, s_j).

### New features / options
| parameter | default | what it does |
|---|---|---|
| `KIARART:photon_escape_fraction_FUV` | required if a FUV group exists | Escape fraction for the FUV group (72e67739d) |
| `KIARART:dust_opacity_FUV_cm2_per_g` | 9.1e4 | Dust opacity per unit dust mass for the FUV group (9925f6e8b) |
| `KIARART:G0_nH_min_cgs` | 1e-3 | G0 from RT is zero below this n_H, so the IGM gets none (31855229b) |
| `KIARART:FUV_LW_photon_fraction` | 0.208 | LW photon fraction of the FUV band (BPASS, Z = 0.014); scales crackle's Draine-based ISRF H2 rate by 1.84 by default, for RT-derived G0 only (455c00997, 12c66e57e) |
| `KIARART:transport_length_floor_softening` | 1 | Floor the RT transport length at the gravitational softening (15d41b9a3) |
| `KIARART:local_cfl_reduced_c` | 1 | Scale c_red by length_courant/length (floor 1e-3) for softening-floored particles; removes G0 knots in sub-softening particles (15d41b9a3) |
| `KIARART:transport_area_ratio_max` | 10 | Cap on sum\|A\|R/V (sphere = 3) in the exchange, time step and c_scale; replaces the compile-time `rt_area_ratio_max` (c11f72e4a) |
| `KIARACooling:G0_computation_method` | — | New value 8: G0 = E_FUV (c_red/c)/u_Habing from the RT FUV group, inside and outside the subgrid ISM. Stops with an error without KIARA RT or a FUV group (31855229b) |

- BPASS tables are read per band name, one per group. A group 0 entirely below 13.6 eV is the FUV group and uses `Table_FUV` (mean energy 9.18 eV, no cross section). `tools/make_kiara_rt_bpass_tables.py` builds the tables from the BPASS SED grid (72e67739d). With three groups the behaviour is unchanged.
- FUV absorption is dust-only. There is a warm-phase path with (1 − f_cold) of the dust, and a cloud path with covering fraction C and chords of 4R/3 (9925f6e8b, 3b37a9109).
- The effective speed of light c_red·max(1, min(ratio, cap)/3) is used for every energy-density↔flux conversion (G0, absorption, photoionization rates); transport, time step and the flux limiter keep c_red. New snapshot field `RTEffectiveLightspeed` (c11f72e4a).
- Outside the ISM, cooling passes −1 − G0 as crackle's `isrf_habing` (decoded by crackle). With other G0 methods this is still exactly −1 (31855229b).

### Behaviour changes
- `KIARART:ism_cold_cloud_temperature_K` is removed (3b37a9109, e9515cc9a). Subgrid clouds no longer go transparent above 300 K. For FUV, only dust absorbs, independent of subgrid temperature; before, FUV heating could make the ISM transparent (G0 ~ 1e7 in dense gas by z = 6). For ionizing groups, clouds absorb at any temperature over 4R/3 chords, up to their case-B recombination rate α_B n_H,cl² φ V dt. The warm phase stays transparent to ionizing photons.
- The LW factor applies only to RT-derived G0 (method 8). Recipe G0 keeps the assumed Draine rate (12c66e57e).
- The KIARA-RT crackle instance now follows the `KIARACooling` dust/H2/ISRF settings (455c00997, see above).

### Other
- Debug outputs (`--enable-rt-debug` / `SWIFT_KIARA_RT_DEBUG`) for checking photon conservation: `RTDebugPendingTransport`, `RTDebugOneSidedTransport`, `RTDebugTransportAccumulated`, `RTDebugEnergyBudget`, `RTDebugEmissionCarry`, `RTDebugEmissionPending`, plus a rate-limited over-export report (15d41b9a3, 08e13c44b, 4b8cfc91c, f76bae273).

### Upstream SWIFT
- FOF speed-ups (b3d5e935d). `FOF:group_id_offset` is removed (an error if set to anything other than 1). `FOF:min_group_size` must be ≥ 2, and the linking-length and `group_id_default` checks are stricter. Restarting at the end of a run now does a full `engine_rebuild` so FOF works over MPI.
- Mesh potential no longer includes each particle's own smoothed mass; the exact-gravity check adds the self periodic-image term (0fd5a1089).
- Stars ghost-loop and self-recursion improvements, with a debug field `N_density_exact_edge` in EAGLE sparts (4d85fd802).
- Small fixes (d659e0a8a): sign of the spart `x_diff` displacement update in `runner_others.c`, default values of optional I/O fields written at the correct type width, GEAR sink fixes, and configure checks for GEAR sink and star-formation dependencies.
- Distributed mesh send-count loops stop at the first empty bucket (4b3b4ddd6). CI python jobs run in a venv (6547ac217).

---

<a id="pr-52"></a>

## PR #52 — Firehose momentum-exchange fixes (2026-11-12, merge 7c533e9e1)

Three fixes to the firehose wind model: the momentum exchange was being thrown away each step, the ambient sums were used unnormalised with the model off, and streams mixed with dense ISM gas. Also makes the MAGMA2 conductivity alpha a runtime parameter.

### Bug fixes
- **Firehose momentum exchange was discarded every step while its dissipated energy was kept** (a937619c4). Since 781844de5 the velocity change went only into `p->v`, which kick2 resets to `xp->v_full`. Streams therefore never slowed down or recoupled on the Mach criterion, and kept reheating their ambient neighbours: an unbounded energy source that pushed CGM gas out of haloes. The change is now applied to `xp->v_full` and the gpart's `v_full` too.
- **With `use_firehose_wind_model` off, the ambient sums grew without bound and were still used** (f88902d7f). The jet recoupling Mach test used the unnormalised `v_ambient`/`u_ambient`, and decoupled winds got a mixing-layer cooling solve every step, skipped dust sputtering, and reported the mixing-layer time as `CoolingTimes`. The sums are now always reset in `chemistry_init_part` and zeroed after the density loop when the model is off. No change with the model on.

### New features / options
| parameter | default | what it does |
|---|---|---|
| `SPH:conductivity_alpha` | 0.075 | MAGMA2 artificial conductivity alpha, previously the compile-time `const_conductivity_alpha`. Printed at start-up and in snapshots (d8a38d963) |

### Behaviour changes
- **No firehose mixing with ambient gas denser than `firehose_nh_ambient_max_cgs`** (9655d03b6). Before, this parameter only capped the density used for the mixing-layer cooling time, not the exchanged mass. Once the momentum fix was in, ISM-launched streams mixed with dense ISM gas, stalled and recoupled inside the galaxy, which roughly doubled the stellar mass. A stream now crosses the ISM as a plain decoupled wind and mixes only with more diffuse gas. The cap on the cooling density is unchanged.

---

<a id="pr-51"></a>

## PR #51 — Size-resolved dust, WN07 SF-law knobs, AGN-wind diagnostics (2026-11-10, merge a76aaad84)

Brings in the grain-size-resolved dust model from PR #49 (Balu Sreedhar): crackle evolves the size distribution, grains collide in the ISM, and dust moves between particles with its sizes. Also adds redshift/density knobs for the Wada-Norman SF law, wind-speed H(z) options, ADAF/quasar-wind options, and a recoupling log.

> ⚠️ Restart files from earlier builds are not compatible.

### Bug fixes
- **SNe deposited on gas were not conserved** (0a7346cc7, 4107a054f). Gas multiplied a summed SN rate by its own step, so SN dust destruction and SN-based G0 were undercounted or overcounted depending on the step ratio. Split particles also each kept all the SNe. Stars now deposit SN numbers, gas uses them at a constant rate over its step (RT sub-cycles included), and splits share them.
- **Dust in gas that skips the cooling solve was never sputtered** (a562b7ee6). Decoupled non-stream gas and AGN-heated gas with cooling shut off are now sputtered in `firehose_cooling_and_dust()`.
- **Carbon-rich AGB test used C > O by mass** (ddce1c474). It is now C > 0.75 O, i.e. C/O > 1 by number. This changes AGB dust with or without bins.
- **BHs did not accrete dust** (9f7072d4b). Nibbling raised the gas D/G and lost the dust's metals. BHs now take the nibbled fraction of the dust, all of it when swallowing, and include it when seeded from gas.
- **MAGMA2 h could drift towards zero and overflow rho to inf** (803931087). The drift exponent is bounded so that h stays at or above `h_min`, and rho is scaled to conserve rho h^3.
- **WN07 epsc redshift scaling leaked into the BH fuel** (e4446f70c). `dense_gas_fraction` now uses epsc without the a^k factor. In a k=1 test this had made BHs 0.5-0.9 dex smaller at z>=2.
- Sputtering time-scale now uses T in K, so it is correct for any unit system; no change for KIARA runs (aacd841e1).

### New features / options
| parameter | default | what it does |
|---|---|---|
| `--with-kiara-dust-size-bins=N` (configure) | 10 | Grain-size bins per dust material (carbonaceous, silicate). 0 builds without size-resolved dust (90e48ef88, c565e2d1d) |
| `KIARACooling:dust_grain_size_min_micron` / `_max_micron` | 0.001 / 1 | Logarithmic grain-size grid |
| `KIARACooling:dust_grain_sizes_micron` | — | Explicit grid of exactly N radii |
| `KIARACooling:dust_initial_source` | SNII | Size spectrum for dust without a size history (SNII or AGB) |
| `KIARAFeedback:dust_{SNII,AGB}_a0_micron`, `_sigma`, `_p` | 0.1; 0.6/0.47; 0/4 | Li+21 stellar dust source spectra |
| `KIARACooling:dust_coagulation` / `dust_shattering` | 1 / 1 | Grain-grain collisions in the ISM (Li+21; Hirashita & Yan 2009) (397c76b0a) |
| `KIARACooling:dust_grain_density` | 2.4 | Grain material density, g cm^-3 |
| `KIARACooling:dust_collision_velocity_cap` | 1 | Caps grain velocities at this many times the gas velocity dispersion across the kernel |
| `KIARAStarFormation:wn07_tff_exponent` | 1 | SFR per dense gas ∝ 1/(t_ff^q t_ref^(1-q)), q in [0, 1] (c836c921f) |
| `KIARAStarFormation:wn07_tff_ref_density_cm3` | 500 | Reference density for t_ref (rho/m_H) |
| `KIARAStarFormation:wn07_epsc_redshift_exponent` | 0 | WadaNorman epsc × a^k before the [0.001, 0.1] clip (571e7c401) |
| `KIARAFeedback:FIRE_velocity_Hz_exponent` | 1/3 | Power of H/H0 in the wind circular velocity (3c2dfba60) |
| `KIARAFeedback:FIRE_velocity_Hz_max_boost` | 0 (no cap) | Caps the extra factor (H/H0)^(exponent − 1/3) |
| `KIARAFeedback:bh_recouple_log` | 0 | Writes a `BH_RECOUPLE` line per recoupling of a BH-launched particle: type, reason, flight time and distance, speed, n_H, T (2426b9b36) |
| `ObsidianAGN:adaf_wind_min_speed_km_s` | 0 | Whole-kernel ADAF: shares too small for this speed are kicked at exactly this speed stochastically, conserving expected energy (d524f1029) |
| `ObsidianAGN:quasar_wind_kicks` | 1 | 0 skips quasar/slim-disk wind kicks but keeps their mass loading in f_accretion (diagnostic) (47f0f0295) |

- **Size-resolved dust** (3074b15be, 5c6cd4ef6, c5f488b68, 0f573e947, 397c76b0a): each gas particle carries a material × size-bin distribution. crackle grows, shocks and sputters grains by size, and scales the surface processes by the grain-area factor. Dust moves with its sizes through stellar injection, diffusion and firehose mixing. Bins are always rescaled to the element dust masses. New snapshot fields `DustCarbonaceousBinMasses` and `DustSilicateBinMasses`. With 10 bins the gas particle is 1184 B (1024 B without bins). `--enable-crackle-validate` is refused when bins are on.

### Behaviour changes
- Size-resolved dust is on by default (10 bins); build with `--with-kiara-dust-size-bins=0` for the previous single-size model (c565e2d1d).
- The AGB C/O fix changes results with default parameters.

### Other
- Dust unit tests (`testKIARADustGrid`, `testKIARADustSources`, `testKIARADustEvolution`, `testKIARADustCooling`) run in `make check`, plus a sanitizer runner `tests/testKIARADustGrid.py`. KIARA dust documentation added under `doc/RTD` (90e48ef88, d99c7dcea).

### Upstream SWIFT
- Analysis scripts updated for the latest swiftsimio; README/script updates; UNIX line endings enforced (22feb14e4, bdf71da05, 172ccbd71).

---

<a id="pr-48"></a>

## PR #48 — KIARA_25 example, unified SF efficiency, optional photometry (2026-11-04, merge 6a85aa342)

Adds a ready-to-run KIARA example in the EAGLE 25 Mpc low-resolution volume, merges the star-formation normalisation parameters into one, and makes the stellar photometry tables optional.

### New features / options
| parameter | default | what it does |
|---|---|---|
| `Stars:luminosity_filename` | none (now optional) | If missing or `none`, no photometry tables are read and the stars' `Luminosities` field is not written (KIARA stars only). If set, nothing changes (15a067b05) |

- **`examples/KIARA/KIARA_25`** (50bd0f5e9): runs the current KIARA calibration (MAGMA2, KIARA cooling/chemistry/stars/feedback, Obsidian BHs) from the DM-only EAGLE L0025N0188 ICs. Their gas particle mass (1.45e7 Msun) is close to the KIARA calibration resolution. `run.sh` fetches the ICs and the Cloudy table; the chem5 tables are in `examples/KIARA/chem5`. `kiara_25.yml` is meant to follow the calibration.

### Behaviour changes
- **`KIARAStarFormation:star_formation_efficiency` now applies to all SF models** (0409ad91f). It is the efficiency per free-fall time for SchmidtLaw and a multiplier of the model SFR for WadaNorman and lognormal. It is now optional with default 1.
- **`KIARAStarFormation:lognormal_epsilon` removed** (0409ad91f). A parameter file that still has it stops at start-up with a message to rename it, so that an old, unused `star_formation_efficiency: 0.5` is not silently picked up. The WadaNorman start-up message now prints the efficiency.

---

<a id="pr-47"></a>

## PR #47 — Built-in crackle, H2/ISRF and dust physics, Sigma_SFR winds, BH dynamical friction (2026-11-04, merge 72aa91149)

A large merge. crackle is now compiled into SWIFT instead of linking Grackle, and runs directly on the particles. It also adds H2 photodissociation/shielding options, a new G0 method, dust-model corrections, chem5-based wind energy, Sigma_SFR-dependent winds, BH dynamical friction, runtime MAGMA2 neighbour numbers, and a large trim of the particle structs.

> ⚠️ Restart files from earlier builds are not compatible.

### Bug fixes
- **Total metallicity dropped chem5's untracked metals** (3f4511ec4, 9e89a7fdc). Z was recomputed as the sum of the tracked elements after feedback, diffusion, mixing, sputtering, jets and SF, so metals in untracked elements (Na, Al, Ar, ...) were lost from Z and from metal cooling. Z is now tracked on its own: feedback adds chem5's total metal ejecta, it diffuses as its own scalar, firehose exchanges it (`dm_Z_total`), and returned dust adds its whole mass. Self-enrichment also keeps these metals in Z.
- **chem5 SN energy counted twice with `wind_energy_from_chem5`** (905fed8d7). Energy fed to the wind reservoir was also deposited thermally when `use_SNII/SNIa_energy_from_chem5` was set. It now goes to the wind while the star is launching, and is deposited thermally only once the wind is complete.
- **SN dust destruction depended on galaxy-scale density** (3375b25a8). The destroyed dust was 22 Msun/Myr times 1/sqrt(G rho), clearing about 40-200 Msun of gas per SN, roughly 10x below standard estimates, and `dust_destruction_eff` cancelled out. Now uses the McKee (1989) Sedov-Taylor shocked mass, 6800 `dust_sne_coeff` (100 km/s / v_s)^2 Msun, with `dust_destruction_eff` again the destroyed fraction.
- **N and Ne dust grew from self-enrichment seeds** (4c9aacadc). Both are dropped from the dust elements; dust arrays now hold C, O, Mg, Si, S, Ca, Fe.
- **Denormal dust mass overflowed to inf** (4c5102e2d). Dust diffusion, firehose mixing and stellar dust injection now test against `FLT_MIN` instead of 0.
- **Crash at step 0 from foreign `stars_sort` tasks sorting uninitialised sparts** (d94c4ed81). They are now skipped in `engine_skip_force_and_kick()`, and stale foreign stars-sort requests are cleared.
- **MAGMA2 neighbour-number check never ran** (f3c1177fb). It was guarded by an undefined macro. Superseded by 8d935ebc6.

### New features / options
| parameter | default | what it does |
|---|---|---|
| `KIARACooling:H2_isrf_dissociation` | 0 | G0-driven Lyman-Werner H2 photodissociation, self- and dust-shielded (dd7723591) |
| `KIARACooling:H2_isrf_dissociation_rate_cgs` | 3.4e-11 | Photodissociation rate per Habing unit, s^-1 |
| `KIARACooling:H2_dust_shielding_cross_section_cgs` | 1.9e-21 | Dust LW shielding cross-section per H at local D/G (Sternberg+14), scaled by the particle's D/G |
| `KIARACooling:H2_shield_mean_density` | 0 | H2 shielding column at the particle's mean, not cold-cloud, density |
| `KIARACooling:H2_sternberg_limit` | 0 | Limits the cold-phase f_H2 to the column behind the Sternberg+14 atomic layer, returning excess H2 to HI (985b3ce88) |
| `KIARACooling:G0_computation_method` = 7 | — | New G0 mode: f_loc × local SFR-density G0 + (1 − f_loc) × a galaxy-scale field from the host SFR, MW-normalised (9c5ed9f01) |
| `KIARACooling:G0_local_fraction` | 0.5 | f_loc for method 7 |
| `KIARACooling:min_subgrid_number_density_cm3` | 0 (off) | Floor on the subgrid cold-cloud density, rho/m_H in cm^-3; seen by every user of that density (826099c1c) |
| `KIARACooling:dust_oxygen_stoichiometric_cap` | 0 | Caps dust O at the oxides of dust Mg, Si, Ca, Fe; excess back to gas (56383a390) |
| `KIARAFeedback:wind_energy_from_chem5` | 0 | Wind energy budget = chem5 SNII (incl. hypernova) energy + 1.3e51 erg per SNIa, without the metallicity energy boost (a7bdae778) |
| `KIARAFeedback:vwind_sigma_sfr_slope` | 0 (off) | Wind speed × clip((Sigma_SFR/pivot)^slope, min, max), Sigma_SFR from the host galaxy's SFR and SF radius (b6cc4391d) |
| `KIARAFeedback:vwind_sigma_sfr_pivot_Msun_yr_kpc2` | 0.1 | Pivot Sigma_SFR |
| `KIARAFeedback:vwind_sigma_sfr_min_factor` / `_max_factor` | 0.5 / 2 | Clip range of the boost |
| `KIARAFeedback:vwind_sigma_sfr_min_radius_kpc` | 0.25 | Floor on the projected SF radius |
| `KIARAFeedback:vwind_sigma_sfr_eta_slope` | 0 | Mass loading × boost^-eta_slope (2 keeps eta v^2 fixed); needs a non-zero slope (7b399532c) |
| `FOF:galaxy_min_sf_particles` | 10 | Galaxies with fewer SF particles get SF radius 0 |
| `ObsidianAGN:dynamical_friction_model` | 0 | 1: Dubois+13 gas dynamical friction on BHs (as in DARWIN-1) where kernel gas is denser than `df_min_nH_cgs`; repositioning below it (caf7d9547) |
| `ObsidianAGN:df_min_nH_cgs` | 0.1 | Density above which friction replaces repositioning |
| `ObsidianAGN:df_boost_exponent` / `df_boost_nH_cgs` / `df_boost_max` | 2 / 0.1 / 0 (no cap) | Boost alpha = max[(rho/rho_crit)^beta, 1] |
| `ObsidianAGN:df_coulomb_log` | 3.2 | Coulomb logarithm |
| `ObsidianAGN:df_use_dynamical_mass` | 0 | Use the dynamical rather than subgrid BH mass |
| `ObsidianAGN:reposition_radius_from_softening` | 1 | Maximum reposition distance = `max_reposition_distance_ratio` × baryon softening (EAGLE); 0 keeps the h_BH-based radius (6ec21a73e) |

- New snapshot fields: `GalaxySFRRadii` (gas, stars), `GalaxyCentreDistances` (gas), `DynamicalFrictionAccelerations` (BHs); named per-element columns for `DustMassFractions`.
- **Built-in crackle** (fd3b14b17, 1a022f808, a92af28d4, 680ee2641, 8bfbb8751, 5f2bc046b): crackle and the Grackle C it needs are imported into `src/cooling/KIARA/crackle`, with symbols in a `crk_` namespace, and built as `libcrackle`. No external Grackle is linked (`--with-grackle` is ignored for KIARA, and KIARA-RT now requires KIARA cooling). crackle now fills its state straight from part/xpart. UVB rates and units are computed once per step, not per particle.
- **`--enable-crackle-validate[=direct]`** (124f71f7a, b7d2a34e7, 544456cc0): bit-for-bit comparison of the built-in crackle against an external Grackle library, or of the direct path against the field-array path.
- **Compile-time element lists** (e917e2403, 40f3ee9cf): `KIARA_ELEMENTS` and `KIARA_DUST_ELEMENTS` in `chemistry_struct.h` set the tracked and dust elements (any of chem5's 32). Initial abundances of unlisted metals default to zero.

### Behaviour changes
- **MAGMA2 neighbour number is runtime-only** (596f4eef9, 8d935ebc6). The default moved from 114 to 64 neighbours in 3D, and then the compile-time `const_kernel_target_neighbours` was removed. `SPH:target_neighbours` (or `SPH:resolution_eta`) is now the only setting. The slope limiter's eta_crit, the low-neighbour fallback threshold and KIARA's self-shielding length floor (0.516 mean interparticle spacings) are derived from it at run time. Parameter files should use `SPH:target_neighbours: 64`.
- **`KIARACooling:dust_growth_tauref` no longer read** (dd963aa82). The growth time at `dust_growth_densref` is fixed at 1 Gyr. A value of 1 is accepted; any other value stops the run and prints the equivalent densref.
- **BH repositioning radius is again in softening units by default** (6ec21a73e). In an s12.5n128 test to z=2 this gave about 5x more repositions.
- Star per-step dust yield is now stored as float, so results are not bitwise identical (a3f1c1dee).

### Other
- Struct trimming, no physics change (a16cb0a11, 3a0c24157, ab643081c, a3f1c1dee): removed write-only cooling/feedback fields and EAGLE tracking fields; the shear tensor and dZ_dt now share a union; stars get a 48-byte birth-composition chemistry struct. KIARA part goes from 1120 to 992 B, xpart from 160 to 128 B, spart from 800 to 416 B.
- `config.h` is included as `<config.h>` so out-of-tree builds work (9556819bf).

### Upstream SWIFT
- Parallel I/O: fixed pointer advance for non-gas fields with conversion functions (e1a963fa6).
- Tree rebuild condition from the stars' neighbour search altered; more compact `pcell_step` layout (7c7a1bacf, 7a1c598ff).
- More accurate SPH kernels (less premature rounding); driven turbulence terms (`TurbulenceDriving:*`) (c56ca8895, 589892dee).
- Cosmology table-difference operators guarded against negative rounding (707af2e81).
- GEAR chemistry: fractional SN counts passed to the yield setter; GasFreeSinkCells example (d93de34e2, 3a17823da).
- Docs: IsolatedGalaxy README, LaTeX onboarding guide (d2dba1ec2, dc9060572, c71451f8a).

---

<a id="pr-46"></a>

## PR #46 — Jet-recoupling runaway heating fix (2026-10-24, merge f44218be8)

Fixes the energy non-conservation that produced T > 1e10 K gas and dt_min crashes in the s25n256 jet/spin calibration runs, plus build fix-ups after the zoom merge.

> ⚠️ Restart files from earlier builds are not compatible.

### Bug fixes
- **MAGMA2 SPH fallback left ill-conditioned pairs without any force** (888025d1c). `G_ij` and `G_rad_ij` started at 0 and were set in only one branch, so those pairs had no pressure, viscosity, PdV or dh/dt. They now start at `dx`. The fallback is also forced for particles at h_max with fewer than half the target neighbour number (new debug field `LowNeighbourFallbackCounts`).
- **MAGMA2 viscosity could reverse a pair's approach and create energy** (888025d1c). Over one step, viscosity may now remove at most half of a pair's radial relative velocity (`const_viscosity_max_dv_fraction`). The viscous acceleration and heating are scaled together, so pairs stay conservative. The step length comes from a new `dt_kick_therm` field set in kick_extra. The cap binds on about 1e-4 of force pairs.
- **Jet particles recoupled inside dense gas at ~5000 km/s** (888025d1c). After its timer expires, a decoupled jet particle now recouples only below a density or Mach threshold (see below). A z=0.5 to 0 restart test had no runaway: final gas thermal energy 2.9e9, against 1.0e10 without the fixes.

### New features / options
| parameter | default | what it does |
|---|---|---|
| `KIARAFeedback:jet_recouple_nH_cgs` | 0.13 | Expired-timer jet particles recouple only below this n_H (or below the Mach threshold) |
| `KIARAFeedback:jet_recouple_mach` | 10 | Expired-timer jet particles recouple only below this Mach number relative to the ambient gas (or below the density threshold) |

### Other
- Removed the unused "top 10 hottest particles" debug block in `engine_step`, which looped over all gas particles every step (888025d1c).
- Zoom-merge fix-ups: gcc 12 `-Werror` fixes in `serial_io.c` and `zoom_init.c`, whitespace after backslashes in `Makefile.am`, and removal of the dead `src/engine_marktasks.c` (fce860ae7).

### Upstream SWIFT
- Sinks: Boss & Bodenheimer and turbulent collapsing-cloud examples; tracer structs added to the other sink modules (411ef6a00, 4b81bf714).
- Hydro vector interaction functions: argument qualifiers updated (7bc846327).
- Clearer description of the progenitor particle IDs (5b87eb066).

---

<a id="pr-44"></a>

## PR #44 — Obsidian BH spin evolution; H2 across subgrid ISM exit (2026-10-21, merge 28d4bdf9a)

Adds black-hole spin evolution to Obsidian (accretion, warps, mergers, jets along the spin axis, spin-dependent efficiencies) with a set of ADAF spin/jet options, keeps H2 consistent when gas leaves or re-enters the subgrid ISM, outputs crackle's own dust temperature, and fixes the post-snapshot `runner_do_stars_sort` segfaults.

> ⚠️ Restart files from earlier builds are not compatible.

### Bug fixes
- **Segfaults in `runner_do_stars_sort()` after a snapshot** (94a8d2c99). Dumping a snapshot, STF output or restart (with `free_foreign_during_restart`) re-links foreign particles and resets each foreign cell's `parts_rebuild`, but other ranks keep sending counts relative to their own last rebuild, so stars formed since then were counted twice and foreign progeny pointed at the wrong sparts. This could also silently pair the wrong foreign stars with local gas. The engine now forces a rebuild right after the re-linking.
- **DustTemperatures were 2.18x too high** (25b38957d). They came from Grackle's `calculate_dust_temperature()`, whose opacity is per gram of gas while crackle's is per gram of dust. Crackle now returns its own Tdust via the new `grackle_field_data` `dust_temperature` field (also in the RT path). Requires crackle b401986.
- **AtomicHydrogenMasses ignored the cold fraction in the subgrid ISM** (25b38957d). In the ISM the species fractions are cold-phase values, so HI mass is now multiplied by `subgrid_fcold`, as for H2; single-phase gas is unchanged.
- **BH spin-evolution bugs** (0945b1ab0, a5f01ff76). The gas angular momentum accumulation was commented out, so the spin direction never followed the accreting gas and the King et al. (2005) test never triggered. The merger remnant spin now uses Barausse & Rezzolla (2009) Eq. 6 (equal-mass non-spinning merger gives a_fin = 0.69, was 0.009). The merger orbital plane uses the periodically wrapped separation. Spin-up uses the Novikov-Thorne efficiency (thin disks reach 0.998, not ~0.97). The thick/slim-disk warp angular momentum was 2x too large. The redirection time-step is guarded against sine = 0.
- **Spin evolved with reservoir inflow instead of mass reaching the BH** (baa1404cc). With the accretion-disk reservoir the spin did not change while the reservoir drained (35-47% of ADAF growth in tests). It is now evolved after the reservoir step, and the reservoir tracks its Bondi-fed mass (`accretion_disk_bondi_mass`) so `bondi_fraction` refers to the draining mass.

### New features / options
| parameter | default | what it does |
|---|---|---|
| `ObsidianAGN:subgrid_seed_spin` | 0.01 | Spin assigned to BHs at seeding |
| `ObsidianAGN:alpha_acc` | 0.2 | Viscous alpha of the subgrid accretion disk |
| `ObsidianAGN:TD_region` | `B` | Shakura & Sunyaev (1973) thin-disk region (`B` or `C`) |
| `ObsidianAGN:fixed_spin` | -1 (now optional) | If set, all BHs use this spin for radiative and jet efficiencies (spin still evolved as a diagnostic); if omitted, efficiencies follow each BH's evolved spin |
| `ObsidianAGN:jet_launch_dir` | 4 (now optional) | New option 4: bipolar jets along the BH spin axis |
| `ObsidianAGN:adaf_torque_phi_factor` | 1 | Magnetic flux of the torque-fed ADAF inflow relative to a MAD; sets its jet efficiency and spin-down (Bondi share is always MAD) |
| `ObsidianAGN:adaf_jet_free_energy` | 0 | If 1, ADAF jets use the full MAD efficiency for all inflow while spin-down still uses `adaf_torque_phi_factor` |
| `ObsidianAGN:adaf_freeze_spin` | 0 | If 1, ADAF accretion does not change the spin magnitude (direction, sign and merger spins still evolve) |
| `ObsidianAGN:adaf_f_accretion_reference_spin` | -1 (off) | If > 0, the ADAF f_accretion uses the MAD jet efficiency at this spin instead of each BH's; 0.9 reproduces the old `fixed_spin` f_acc |
| `ObsidianAGN:torque_accretion_max_radius_softening` | 0 (whole kernel) | Only cold gas within this many Plummer-equivalent baryon softenings feeds torque accretion and sets the disk direction |

- BH spin evolution (nicolethom: 0e7dc58d4 ... 92d6c09d1; Romeel Dave: bb56a0bb2, a349bd8c1 and above): new `black_holes_spin.h` with thin/slim/thick-disk spin-up, disk warps, King et al. alignment, merger remnant spins, slim-disk jet spin-down (Ricarte+23) and Narayan+22 MAD spin-down for the Bondi share of ADAF inflow; torque-fed ADAF share treated as a thick disk (0.45 L_ISCO). Spin integrated in sub-steps (|da| <= 0.01, dlnM <= 0.05); a retrograde spin spun through 0 flips. Radiative efficiency (Madau+14/Lupi slim-disk fit) and MAD jet efficiency (Tchekhovskoy+10 with Narayan+22 phi(a)) are per-BH. `adaf_f_quasar_psi > 0` requires `fixed_spin`.
- Disk direction follows `torque_accretion_method` gas selection (`angular_momentum_accreting_gas`).
- New BH outputs: `Spins`, `AngularMomentumDirections` (also optional ICs inputs), `RadiativeEfficiencies`, `JetEfficiencies`; spin, jet efficiency and `fbondi=` in BH_DETAILS; spins in BH_MERGER.

### Behaviour changes
- `jet_launch_dir` default is 4 (spin axis) when omitted; `fixed_spin` is no longer required, and omitting it makes efficiencies spin-dependent.
- H2 kept across subgrid ISM exit/re-entry (1e1c3b8f3, 0dbcde685): on exit (cooling, wind launch, BH feedback, full evaporation) only the cold phase's share of H2 and H2+ is kept, the rest returning to HI (+HII for H2+); on entry the carried H2 is assigned to the cold phase (`H2_cold = min(H2/fcold, X_H)`) instead of reset to 0. MolecularHydrogenMasses outside the ISM is now `m x H2` (single-phase). Requires crackle "Evolve H2 outside the ISM".
- `cooling_subgrid_exit()` now performs all ISM-exit resets (`subgrid_temp`, `subgrid_dens` to the physical density, `subgrid_fcold`, `G0`) at all nine exit points; `feedback_recouple_set_flags` takes the xpart (0dbcde685).

### Other
- Removed unused Obsidian disk properties (`mdot_crit_ADAF`, `edd_crit_thick`, `cs_0_ADAF`, `h_0_ADAF`, `beta_acc`, `gamma_acc`, etc.) and `delta_ADAF`; `spec_angular_momentum_gas` replaced (a5f01ff76).

### Upstream SWIFT
- Limit MPI tags to the MPI standard upper bound (ca6ea3d4b).
- Fix gpart bug when a star is inhibited (3be300fb4).
- Documentation of the timeline and time-stepping scheme with figure scripts (bf00bb5ee).

---

<a id="pr-43"></a>

## PR #43 — RT sub-cycling deadlock and photon-conserving transport (2026-10-18, merge f5aebdc4d)

Three KIARA RT fixes: an MPI deadlock in RT sub-cycling, asymmetric pair exchange in the RT transport sub-cell loops, and a positivity limiter that no longer creates photons.

> ⚠️ Restart files from earlier builds are not compatible.

### Bug fixes
- **MPI deadlock in RT sub-cycles from stale RT times on foreign cells** (9b3ba420e). Foreign cells' RT times are only updated from the super level down, so a foreign top-level cell could look RT-inactive while some of its progeny were active. `engine_unskip` and `engine_unskip_rt_sub_cycle` pruned on the top cell, so a send was never activated while the other rank posted the matching recv. This hung s25n256 on 32 ranks around step 1200. The unskip now walks below foreign cells above their super level regardless of RT time (`engine_rt_cell_time_is_stale`). The run that hung at steps 1249/1302/1319 passed step 1890 with the fix.
- **RT transport lost one side of some pairs in the sub-cell recursion** (b42251829). With depth-limited recursion, a pair whose particles live at different levels was exchanged in two non-symmetric calls, and the smaller particle's side was lost when the larger one was not in an adjacent cell (~1.6% of such calls; 5-10% of injected energy created by z=8.5). For `TASK_LOOP_RT_TRANSPORT`, `DOSUB_PAIR2` and `DOSUB_SELF2` now interact all particles at the level where some particle no longer fits the daughter cells and stop recursing there. Pair-level energy exchange is exact to round-off, and a single-rank s12.5n128 run to z=8.5 is 8% faster. The hydro force loop has the same structure (~0.06% momentum loss in such pairs, MAGMA2, z=9) and is left unchanged.
- **RT positivity limiter could still create photons** (efd0eb719). The per-interface cap used stale or incomplete information, so particles exported more than they held and the clamp to zero created energy (up to 4.5x the injected energy in the c05 tests). The cap is now `E V (A_ij / sum A) dt_frac`, where `sum A` is the larger of the last transport step's and this step's gradient-loop interface areas (new `rt_interface_area` helper and `rt_part_data.area_sum_grad` field), `dt_frac = 2^(bin_min - bin_exp)` limits exporters with longer steps, and a particle with no transport-step area yet (e.g. a recoupled wind particle) exports nothing. Clamped energy drops from 0.25-0.74 to <= 0.003 of the injected energy (s12.5n128, `f_reduce_c` 0.05, to z=5).

---

<a id="pr-42"></a>

## PR #42 — KIARA RT transport overhaul and variable reduced speed of light (2026-10-16, merge c2215d710)

Merges the `rt_fix` RT development line and a series of fixes that make KIARA RT conserve photon energy (correct comoving light speed, BPASS metallicities, volume-change and drift corrections, positivity) plus a variable reduced speed of light and a cold-cloud absorption model for the subgrid ISM. Also carries the wind energy-reservoir deficit, a stars-sort check, and debug instrumentation from the `rt_fix` branch.

> ⚠️ Restart files from earlier builds are not compatible.

### Bug fixes
- **Radiation moved a factor a too slowly in comoving coordinates** (0ce1e17a2). The transport updated comoving quantities over a physical step, i.e. light moved at c_red comoving (7-10x too slow at z=6-9). The update is now scaled by 1/a, and the Courant step becomes `a * psize / c_red` with psize the particle volume radius (as GEAR), which is stable at all redshifts and allows ~3-4x longer RT steps at z >= 6.
- **Wrong BPASS metallicity grid** (6d30ac9bb). From the third entry the grid read 0.01, 0.02, ..., 0.4 instead of 0.001, 0.002, ..., 0.040, so stars got the photon yields of a population up to 10x more metal-poor.
- **Radiation exchanged with only part of the kernel** (685c81e84). `rt_check_ineligible` dropped pairs with r > h although the geometry uses the full kernel (kernel_gamma h), removing ~54% of the transport moment: radiation propagated at ~0.46-0.49 c_red (beam test: 7300 km/s for c_red = 15000 km/s; 15100 km/s after the fix). Pairs with decoupled or cooling-shut-off particles remain excluded.
- **Radiation energy created by volume changes, drift and clamping** (7f5e08273, 251f6e449, dcb9e5f3a). The h-based `(h_old/h)^3` rescaling in `rt_init_part` is replaced by `rt_update_volume`, which rescales energy density and flux by `V_old/V_new` right after the new volume is computed in the ghost (no-op for other RT schemes). The non-conservative radiation drift correction in `rt_predict_extra` is removed (it created ~700x the emitted energy by z=8.5 in s12.5n128). RT rates are now zeroed only in the subgrid ISM rather than for all gas with n_H > 0.13, and each face energy density is capped at twice the particle's own; the old behaviour created ~1000x the emitted energy by z=9.
- **Unstable RT time-step and runaway sub-cycling in dense gas** (6edc9ef9f, 73e83c49e, 4d75f58cf). The RT step now uses the particle's accumulated interface areas, `c_red dt sum|A_ij| / (a V) <= CFL`, with 1/10 of the volume-radius crossing time before the first transport step. The area ratio used in the step is capped at 10 and a conservative flux limiter (each interface exports at most `|A_ij| / sum|A|` of the held energy) keeps those particles positive. The transport length is floored at the baryonic softening (s12.5n128 went from ~480 to ~12 RT sub-cycles per hydro step).
- **Wind kicks could overdraw the SN energy reservoir without charge** (f57ea35ec). The reservoir was reset to zero after the last whole-particle kick overdrew it. The deficit is now kept, and the star launches nothing until later SN energy pays it back.
- **Grackle set u directly for all runs after the rt_fix merge** (ee7a77376). `cooling_do_grackle_cooling` gains a `set_u_directly` argument: the RT thermochemistry passes 1 (set u, zero du/dt), `cooling_cool_part` passes 0 (du/dt from the Grackle energy change, as before). The dt_min error message now lists the time-step limits in label order.
- **Stars-sort merge could segfault on missing or stale progeny sorts** (b5e192d00). `runner_do_stars_sort` now checks that every progeny with stars has an allocated, flagged sort and that counts add up, and otherwise stops with an error giving cell, counts, flags, step and time.

### New features / options
| parameter | default | what it does |
|---|---|---|
| `KIARART:f_reduce_c_max` | `f_reduce_c` | Upper limit of a spatially varying reduced speed of light; equal to `f_reduce_c` means uniform c |
| `KIARART:varc_h_ref` | required if `f_reduce_c_max > f_reduce_c` | Comoving smoothing length in `c = c_min (h / varc_h_ref)^varc_alpha` |
| `KIARART:varc_alpha` | 1 | Power of `h / varc_h_ref` (1 gives equal RT time steps for all particles) |
| `KIARART:ism_cold_cloud_temperature_K` | 300 | Subgrid ISM colder than this absorbs radiation through the cold-cloud covering fraction; warmer subgrid ISM is transparent |
| `KIARART:ism_clouds_per_particle` | 1 | Number of cold clouds N per particle; covering fraction `C = N^(1/3) phi^(2/3)` |
| `KIARART:debug_grackle_heat` | 0 | (only with `--enable-rt-debug`) Second Grackle solve without RT heating to compare delivered vs intended heat |
| `KIARACooling:with_radiative_transfer` | 0 | Passes RT rates to Grackle (`use_radiative_transfer`); forces `with_UV_background` to 0 with a warning |

- Variable reduced speed of light (b61f3b29b): the transport now evolves W = (c E, F) so radiation is not reflected where c changes; the GLF solver works at unit signal speed and each particle scales its F flux by its own c. Uniform c reproduces the old scheme to 4-5 digits. The HLL solver is not converted and now errors. New snapshot field `RTReducedLightspeed`.
- Subgrid ISM absorbs only through the covering fraction of cold clouds, over path `L = c_red dt`: absorbed fraction `1 - (1 - C)^(L / 2 R_p)` (7095d6c70).
- New configure option `--enable-rt-debug` (`SWIFT_KIARA_RT_DEBUG`) compiles the per-particle photon energy budget (`RTDebugEnergyBudget`), `RTDebugInterfaceAreaRatio` and `RTDebugGrackleHeat` outputs (7104cfe3c, 43ab9799f, ae4d8b8df, ce7d951ec, 18f6c00a8, 1f4cc528a, 3b724432b).

### Behaviour changes
- KIARA RT results change substantially with default parameters (light speed, BPASS grid, kernel neighbours, time-step, absorption in the ISM, no drift correction).
- KIARA winds: the SN energy reservoir can now go negative, delaying further launches (f57ea35ec).
- `DOPAIR1_BRANCH_STARS` now warns instead of erroring on interacting unsorted hydro cells (from the `rt_fix` line).
- The per-step `[rt-sc]` RT sub-cycle printout is disabled (the guard was changed to `nodeID == -1`).

### Other
- `rt_fix` branch commits (Zhen Xiang, boson112358; 88d8844d1, 6b1771f47, 0f6738901, 32f03db7f, 387f5e7e7, 703a6a419) and 581b728de, 8946c00f6 add debug instrumentation, still present after the merge. MAGMA2 `struct part` gains ~11 kick/drift internal-energy debug fields, `cooling_part_data` gains ~11 debug u/du fields (wind, SN, cooling, chemistry) and `ncool_this_step` (new output `NumberOfTimesCoolingThisStep`), and `engine_step` runs a "top 10 hottest particles" scan of all local gas every step (the print is commented out). `cooling_copy_to_grackle` also overwrites the x-velocity slot passed to Grackle with the particle ID. The RT timestep now also takes the gravity properties (for the softening floor).
  - *Status today:* the "hottest particles" scan was removed in 888025d1c (PR #46). The MAGMA2 debug fields are now behind `MAGMA2_DEBUG_CHECKS`, which is off by default, and the cooling debug fields are gone. The x-velocity slot still carries the particle ID; crackle uses it only as an ID. The stars loop still warns, rather than erroring, on unsorted hydro cells.

### Upstream SWIFT
- Activate hydro drifts with the rt_gradient loop (575193685, 260272c1b).
- Fix missing vertex weight for pairs with a foreign ci in partitioning (7b65164c2).
- Fix typo in `task_overlap` for sink tasks (b538a89e9).
- AGORA feedback reads the old time bin for the feedback time (a73676a6a); GEAR stellar-wind lookup uses absolute Z, updated feedback table names (4d7a5e2ba).
- REMIX base hydro example changes (5f1eb10c3); FLAMINGO/IdealisedCluster example fixes, EAGLE-XL subgrid (be281ba2a, d97d0280e).

---

<a id="pr-41"></a>

## PR #41 — KIARA review fixes, formatting and documentation (2026-10-12, merge e6e4e63bd)

A small clean-up merge: three correctness fixes in KIARA cooling, feedback and RT, plus clang-format and complete doxygen comments across the KIARA, Obsidian and MAGMA2 modules.

### Bug fixes
- **G0 methods 4 and 5 used one step's SN count instead of the SN rate** (f2ab216cd). `SNe_ThisTimeStep` already holds a rate, but G0 multiplied it by the time-step. G0 now uses the rate directly, with `G0_factorSNe` normalised to the MW rate of 0.015 SN/yr in internal units (as `G0_factor2` is for the sSFR). `cooling_compute_G0` no longer takes the time-step.
- **Oxide oxygen used AGB condensation efficiencies for S, Ca and Fe** (fb7d6e52d). Oxygen locked in oxides took the selected (SNII or AGB) efficiencies for Mg and Si but always the AGB ones for S, Ca and Fe. All five now use the selected table.
- **KIARA RT negative-photon check never fired** (4dd5fe7a0). The BPASS table check tested `fabs(x) < 0`, which is never true; it now works. The 2D stellar-emission interpolation called `exit(1)` (or returned 0 for a NULL table), which could leave other MPI ranks hanging; it now calls `error()` with the same diagnostics.

### Behaviour changes
- `KIARAFeedback:dust_boost_factor` default changed from 2 (Simba value) to 1 (fb7d6e52d).
- With `KIARACooling:G0_computation_method` 4 or 5, G0 from SNe now scales with the SN rate rather than the per-step SN count (f2ab216cd).

### Other
- clang-format 20.1.8 (repository style) applied to chemistry, cooling, feedback, rt, star_formation and stars/KIARA, and to Obsidian, MAGMA2 and the `WITH_FOF_GALAXIES` sections of `fof.c`; whitespace only (2c40e4f04, 445819704).
- Complete doxygen documentation (`@brief`, `@param`, `@return`, struct-member comments) for every function and struct in the KIARA modules, with stale or wrong comments corrected; comments only (8d98ce0dc).

---

<a id="pr-40"></a>

## PR #40 — Physics audit: cooling, firehose, enrichment, star formation and Obsidian fixes (2026-10-12, merge dc942619a)

A large set of correctness fixes across KIARA cooling, firehose winds, chemistry diffusion, stellar feedback, star formation, MAGMA2 handling of decoupled winds and Obsidian BHs. It also adds dynamical-friction BH mergers and an ADAF energy ramp. Many fixes change results with default parameters.

> ⚠️ Restart files from earlier builds are not compatible. (New bpart fields `jet_power` and `df_sat_mass`; new `chemistry_part_data` fields `dZ_dust_dt` and `v_ambient[3]`.)

### Bug fixes
- **Cooling: wrong mu and species conventions** (b1b811d24, 2a34765f0, 3b4d41110, 8015de29b). `cooling_convert_temp_to_u` divided by 1/mu, so a T→u→T round trip inflated T by mu² (~1.5 neutral, ~5 molecular). The electron fraction's He terms were 4x too large. The Ca density pointed at the S slot. Si solar abundance was 10x too high. The warm subgrid phase is now assumed neutral atomic (it no longer takes the cold phase's molecular mu). Subgrid pressure balance now includes both phases' mu (up to ~1.9x higher density for fully molecular cold gas).
- **Cooling: subgrid density bugs** (2a34765f0, f6c988b67). `subgrid_threshold_n_H_inv` was applied twice, so the density-dependent cold fraction never activated. The subgrid density is now the cold-cloud density rho u / u_c, used directly by grackle, RT, SF, G0 and birth density. Before, grackle evolved the cold phase at ~(1−f) of it. The grackle density rises 2–10x. SF/G0 densities are unchanged at f = 0.5 and up to 9x higher at f = 0.9.
- **Dust counted as hydrogen** (96a88a1c5). X was recomputed as 1 − Y − Z_gas in five places, ignoring dust. Dust is now subtracted (X + Y + Z + D = 1).
- **Dust sputtering and self-enrichment** (2a34765f0). Sputtering returned M_d(1−r²) instead of M_d(1−r). Self-enrichment dust fractions now sum to 1.
- **Non-finite checks removed by `-ffast-math`** (b32c18f5c, 548624a92). `fpclassify`/`isnan` checks are now exponent-bit tests. Runs that were silently passing NaN/inf to grackle will now stop with an error.
- **Grackle field struct and fixed redshift** (d9db6b388). The field struct is zeroed, `provide_volumetric_heating_rates` defaults to 0, and a user-fixed redshift is no longer overwritten. New-particle cooling fields are initialised.
- **Metal diffusion off by a⁵ and not conservative** (1b02158b7, aecd04962). The kernel term was comoving: ~250x too slow at z=2. The pair term is now symmetric. The fractional-change limiter, which destroyed metals, is replaced by a cap on D. Dust diffuses per element. `max_fractional_Z_transfer` is still read but unused.
- **`local_sfr_density` was comoving** (e72f5994b). It was used as physical in G0 methods 1, 3 and −3, making G0 too low by (1+z)³. It is now physical and includes the particle's own SFR (521855ef0).
- **Firehose mass exchange used Δln(a) as time** (90a40ba0a). Exchange compared Δln(a) with physical timescales, so the growth term usually won and dm was zeroed. It now uses physical time over the longer step of the pair.
- **Firehose exchange physics** (6078d821b, c0ecbcf28, 4140868ad, 8138dbf97). Fixes:
  - Dust was never exchanged.
  - Dissipated KE used delta_m instead of the reduced mass.
  - The recoupling Mach used the step's dv instead of |v − v_ambient|.
  - Evaporation accounting was wrong.

  Mixing is now capped at `max_fmix` of each particle's mass per step on both sides. This replaces the frame-dependent KE limiter. Ambient gas is also entrained in the growth regime (dm > 0), so streams decelerate.
- **Firehose recoupling and time-steps** (23353397a, e55a5e3b2, e55a634e2, 14d3deeb3). Recoupling is now tested every step, not only after an exchange. Decoupled winds get a Courant limit 2 CFL h / (|v − v_amb| + c_s). The `firehose_nh_ambient_max_cgs` cap now applies to the stream density in the mixing layer too. Expiry of the cooling shutoff no longer recouples winds, which had cut Obsidian BH winds short.
- **MAGMA2 coupled to decoupled winds** (9b8e1cafc, d1f27ab1f). Coupled gas no longer includes decoupled neighbours in its velocity/u gradients and limiter bounds. Approaching winds no longer cut its time-step to h / v_wind.
- **Stellar ejecta over the wrong interval and mass budget** (666770a2e, 57b9dbeac, 574569a8d, 85919d30f). Ejecta were computed for the previous step: overlapping intervals re-ejected mass. With that fixed, the 0.5 ejecta fudge factor is removed (an SSP now returns 0.40 by 13 Gyr). There is no yield extrapolation above 40 Msun and no ejecta above M_u2 (0.390 by 13 Gyr). Net AGB yields may be negative. Negative totals are clamped quietly and conservatively (abf72553a).
- **SNII counted twice in N_SNe** (29e576ce4). The count was 0.021 instead of 0.0103 per Msun for Chabrier. SN wind energy and the grackle SN dust rate halve.
- **SN rate on gas accumulated forever** (5e39fd1c5, 3b4340419). `SNe_ThisTimeStep` is now decayed after each cooling step, including KIARA_RT runs.
- **Feedback mass conservation and kicks** (548624a92, 82cd8fdba, 6a2049cb6, 6012b5115, 4a2bbecc0, 9296d1112). Fixes:
  - Ejecta capped by `max_mass_increase_factor` now stay in the star instead of being lost.
  - Kick directions are seeded per gas particle (all kicks from a star went the same way), and the periodic direction is fixed.
  - `hot_wind_internal_energy` was u(1 K), so hot winds never launched (non-firehose runs only).
  - The stream-radius observational cap was compared in mismatched units.
  - The metallicity eta boost now applies to the wind mass.
  - SNIa evaporation energy accounting is fixed.
- **Star formation** (3fb9bb7e4, 34f42c5b9, 1839f867a, b2f9f7f74, 812ef58e3, 85096e57d, 12e7d3c87, af7d91b0e). Fixes:
  - WN07 used a rate 3.4x too low per t_ff (retune `lognormal_epsilon` 0.2 → ~0.06).
  - The lognormal width and t_ff density now follow WN07.
  - The KMT clumping factor was 30 × parameter (900 by default).
  - KMT f_H2 is set to 0 for s ≥ 2 (it was negative before).
  - The last-SF record is kept when eligible gas has zero SFR.
  - New-star displacements were biased toward −x.
  - New stars now include dust-locked metals.
  - The star time-step floor no longer overrides the maximal step.
- **Obsidian accretion, feedback and RNG** (572f2483c). Fixes:
  - Torque accretion methods 0–2 always gave zero.
  - The ADAF redshift scaling now also applies to the jet efficiency and the loadings.
  - Wind speeds had a double km/s conversion.
  - Dynamical-time method 1 used the wrong a exponent.
  - Correlated random numbers are now independent.
  - The Eddington cap is always applied.
  - Mergers now carry the reservoirs.
  - Jets are kicked at v_jet only, with no jet below the ADAF mass limit.
  - The ADAF radiative efficiency is now continuous with quasar mode.
  - Stellar eta suppression used the velocity suppression redshift instead of `wind_eta_suppression_redshift`.
- **FoF galaxy data** (7b8e51da4). Stale galaxy data on ungrouped particles and a stale `specific_sfr` are now zeroed.
- **Output and debug builds** (f31ac2b43, 2a34765f0, 41d82f650, d9db6b388). BH `MetalMasses` wrote garbage components. `HeIMasses` wrote HeII. Several unit flags and descriptions were wrong. The `FIREHOSE_DEBUG_CHECKS`, `OBSIDIAN_DEBUG_CHECKS` and mode-1 builds now compile.

### New features / options
| parameter | default | what it does |
|---|---|---|
| `ObsidianAGN:bh_df_merger_time_factor` | 1 | Factor f on the Chandrasekhar dynamical-friction time for a secondary BH sharing a FoF galaxy with a more massive BH. After t_df the secondary is moved toward the primary and merged by normal swallowing. ≤ 0 disables (3fb2c3c3e). |
| `ObsidianAGN:bh_df_sink_fraction` | 1 | Fraction by which separation and relative velocity shrink each time t_df elapses; the clock then restarts. ≤ 0 disables (0d96c4f48). |
| `ObsidianAGN:adaf_energy_ramp_width` | 0 (no ramp) | Multiplies the ADAF heating/kick energy by min((M_BH/M_limit − 1)/width, 1). 1 matches the jet ramp (239af5df4). |
| `ObsidianAGN:adaf_max_coupling` | `adaf_coupling` | Upper limit on adaf_coupling·(1+z)^adaf_z_scaling (572f2483c). |
| `SPH:initial_temperature` = −1 | — | New value: start gas at T_CMB(z_begin) (cosmological runs only) (816baadfd). |
| `FOF:cold_gas_temperature_threshold`, `FOF:cold_gas_n_H_threshold_cgs` | 1e5 K, 0.13 | Now used: galaxy gas mass is cold dense ISM gas, excluding winds, instead of star-forming gas (7b8e51da4). |

- BH-driven winds now get a firehose stream radius, so mixing and entrainment apply to them (521ad5470).
- Obsidian reports the kinetic jet power (new bpart field `jet_power`).

### Behaviour changes
- Nearly all fixes above change results with default parameters. The largest effects come from subgrid/grackle densities, the halved N_SNe, the removed 0.5 ejecta factor, the 3.4x WN07 rate, the KMT clumping (900 → 30), physical metal diffusion and G0, and working firehose exchange/recoupling.
- `KIARAFeedback:SNe_smoothing_time_in_Myr` is now used (read by cooling). The default 0 drops the rate after one step.
- `KIARAFeedback:minimum_galaxy_stellar_mass_Msun` is read by the BH model regardless of `suppress_growth`.
- Removed: `KIARAFeedback:eta_suppression_factor_floor` and `KIARAFeedback:galaxy_particle_resolution_count` (never used) (5aeab63e4).

### Other
- Performance: the run-time grackle cooling time is computed only for firehose streams. `CoolingTimes` is computed at snapshot time (88a5e8379). Ambient sums are cheaper, and D = 0 in the IGM skips diffusion (bb79afff5). The Smagorinsky |v| cap is dropped (521855ef0). Heap allocation is removed from the grackle driver.
- An unused KMT calculation is added to cooling (f9b5ee6e5). A call signature is synced with the main repo (1082dd638). The stale `black_holes.h.rad` copy is removed.

### Upstream SWIFT
- Tracer properties are updated when splitting particles (a7fe7ab7b). The gravity epsilon check and its docs are updated (db37e9806).

---

<a id="pr-38"></a>

## PR #38 — Upstream sync with Obsidian accretion/jet rework, FIRE G0 fit and FOFGalaxyIDs (2026-10-08, merge af317b618)

Mostly a sync with upstream SWIFT master (TREECOOL cooling, GEAR BH model, FOF overflow and MPI fixes). The KIARA changes are mainly in Romeel Dave's merge commits 0d2fc2f00 ("Fix merge conflict … GEAR BH model") and 980c791c3 ("… TREECOOL cooling module"), not in standalone commits: Obsidian accretion/jet logic, the KIARA G0 and dense-gas-fraction handling, and dust bookkeeping. The merge message says the PR includes a segfault fix but does not name the commit.

> ⚠️ Restart files from earlier builds are not compatible. `bpart` loses `gravitational_ngb_mass` and `num_gravitational_ngbs`, and the KIARA `cooling_part_data` loses `SNe_ThisTimeStep` (0d2fc2f00).

### Bug fixes
- **ADAF coupling redshift scaling was never applied** (980c791c3). The `min(pow(1+z, adaf_z_scaling), 1)` result was thrown away. The scaling now applies when `ObsidianAGN:adaf_z_scaling` > 0 and is capped at 0.2.
- **Dust mass was not divided when a gas particle was split** (0d2fc2f00). `cooling_split_part` now divides `dust_mass` by n, like `radiated_energy`.
- **Dust condensation did not reduce the star's total ejecta metal mass** (0d2fc2f00). `feedback_dust_production_condensation` now works on `sp->feedback_data.metal_mass` directly, subtracts the condensed dust from `total_metal_mass`, and stores the total new dust in `delta_dust_mass[0]`.
- **Fatal errors turned into warnings in KIARA feedback** (0d2fc2f00). An invalid enrichment weight `Omega_frac` now gives a warning and no feedback to that particle, instead of an error. Negative stellar ejecta and negative ejecta metal masses are now warned about and set to 0. "DUST EXCEEDS MASS" is now a warning, triggered at a dust mass above 0.2 m_gas, instead of an error at a dust mass above m_gas.
- **`FOF:seed_black_holes_enabled` was required even in builds without black holes** (ee1b064f8). It is now forced to 0 under `BLACK_HOLES_NONE`.

### New features / options
| parameter | default | what it does |
|---|---|---|
| `ObsidianAGN:lum_thresh_always_jet_1e45_erg_s` | 0 (off) | Existing option with new meaning. 0 now means off (threshold FLT_MAX). A negative value gives a threshold of \|value\|·10^max(2−z,0), so the threshold rises at low z. Above the threshold a slim-disk BH also launches jets (980c791c3). |
| `ObsidianAGN:dynamical_time_factor` | 0 (off) | Default changed from 1. Stored internally as its inverse. The accretion-disk reservoir (mass fed to the BH over this many dynamical times) is now used only when the value is > 0; before, it was always on (0d2fc2f00). |

- **FOFGalaxyIDs**: builds with `WITH_FOF_GALAXIES` write the snapshot field `FOFGalaxyIDs` ("ID of the baryonic galaxy") instead of `FOFGroupIDs`, for gas, gpart, sink and star particles (46283e8d4).
- **New FIRE-based G0 fit**: `G0_computation_method` 6 uses the D. Salim 5 Aug 2026 fit. It depends on t_ff from the subgrid density, T_warm and the electron fraction, and blends linearly between low-z and high-z forms for 4<z<5. G0 methods 4 and 5 (SN-rate based) are enabled again and read `feedback_data.SNe_ThisTimeStep` (0d2fc2f00).

### Behaviour changes
- **Obsidian gas eligibility** (980c791c3). A new `black_hole_gas_is_eligible()` is used in both swallow and feedback. It skips decoupled or cooling-shutoff gas, BHs with no galaxy or no neighbour gas, and gas in a different galaxy. The boundness cut (u_kin > u_pot) is gone, and feedback now also requires the gas to be in the same galaxy.
- **Obsidian accretion** (0d2fc2f00, 980c791c3):
  - `ngb_mass`/`num_ngbs` now count every neighbour in the density loop, replacing the separate gravitational counters.
  - Torque method 3 now tests `cold_gas_mass > 0` instead of `galaxy_mgas > 0`.
  - Hot/cold gas classification ignores decoupled and cooling-shutoff gas.
  - The ADAF `f_accretion` is now recomputed each call as 1/(1+2ε_jet(c/v_jet)²), using the jet subgrid velocity. The local a-rescaling for negative `jet_subgrid_velocity_km_s` is computed but not used.
- **Obsidian jets and kicks** (0d2fc2f00, 980c791c3). Jet velocity is now added to `v_kick` instead of replacing it. `ObsidianAGN:quasar_lum_thresh_1e45_erg_s` (outward quasar blowout above a luminosity) has been **removed**. Gas kicked or heated by AGN now has `dense_gas_fraction` set to 0.
- **KIARA star formation** (0d2fc2f00). The KMT H2 calculation moved inside the Schmidt-law SFR and no longer requires s<2. KMT now raises an error with any other SF model. The separate Grackle-H2 `dense_gas_fraction` branch is gone, so Schmidt-law SF uses 1 for non-KMT models. WN07 stores `dense_gas_fraction = min(f_c·ε_c/0.01, 1)`.
- **KIARA chemistry**: the diffusion coefficient is set to 1e-10 in gas below the mean baryon density (980c791c3).
- **BH_DETAILS output**: added `MBHold` and removed `mngb`. The torque rate is now printed multiplied by `f_accretion` (0d2fc2f00).

### Other
- Debug output: more detail in the low chemistry-timestep message (later commented out). New MAGMA2 "LOW COURANT/KICK TIMESTEP" warnings when dt_min < 1e-10. Obsidian BH_KICK/BH_HEAT diagnostics moved under `OBSIDIAN_DEBUG_CHECKS` and hard-wired to one BH id (0d2fc2f00, 980c791c3).
- Merge leftovers on master:
  - Unresolved `<<<<<<< HEAD` / `>>>>>>>` conflict markers remain in the `--with-subgrid` help string in `configure.ac`. *Fixed in b897bf792.*
  - `part.h`, `hydro.h`, `hydro_io.h`, `hydro_parameters.h` and `debug.c` gain `MAGMA_SPH` branches for the `src/hydro/MAGMA` scheme. `configure` has no option that selects it.
  - Upstream conflict resolutions in `space.h`, `fof.c` and `configure.ac` (825104867, 7bb15e5a0, b821ca74f).

### Upstream SWIFT
- **New physics/models**:
  - TREECOOL primordial cooling, Katz+96 (3de99061f), with parameters `TREECOOLCooling:TREECOOL_file`, `rapid_cooling_threshold`, `UV_background_start_redshift`, `with_Compton_cooling`, `log10_T_min/max`.
  - GEAR black hole module and `--with-subgrid=GEAR-BH` (9fda2098e).
  - AGAMA external potential, `--with-agama` (7d69828db).
  - AGORA CMB floor, `GrackleCooling:agora_cmb_temperature_floor` (20deec551).
  - GEAR feedback tracers and other GEAR fixes (f5d714cd6, 2ab96996d, 821767233, ccef0f8a1).
- **Correctness/MPI**:
  - FOF overflow fix with custom MPI ops (f4545acaf).
  - Particle-buffer alignment in `proxy.c` (fec059aa7).
  - Hydro sort-array races and limiter locks (aa9c871db).
  - Recompute h_max after limiter unpack (3f08ddded).
  - Star-resort activation depth (dba956ba6).
  - Lock release in stars subset pair (9aa2ff884).
  - Non-periodic gravity interaction counts and faster pair walks (3fb916aee, 098cc9697, d982a70af).
  - Negative time steps removed (f36a024ac).
  - `MPI_BYTE`/`MPI_BOR` in place of `MPI_CHAR`/`MPI_MAX` (e8a174bbf).
  - Spin-hint locks (6f96466a1).
- **Sinks/SF (GEAR)**: e.g. 9e8b7a48c, 124a3cea5, a58ba5000, 50d567668, e730d13e1, 15acd82e7, 6e06c5605.
- **Input validation**: early checks on parameter, file-name and buffer lengths, lightcone/Ewald file errors, CFL/eta/EAGLE parameter checks, and comm-buffer allocation checks (4707b95b0, de319c5ca, 1f40efe76, 449ebbf09, d347cf647, 8221e7469, 931fd6e5f, and others).
- **Docs/tests/build/examples**:
  - `parameter_example.yml` completed (446ab121d).
  - Out-of-source test fixes (ba4b1ff05, aa59c9d12, and others).
  - HBT-HERONS example and DMO-flag removal (f66067764, 9cc766972).
  - Nbit filter for long longs (3aef509bd).
  - AREPO metal-advection and Balsara–Kim examples (e53e77896, 5e5bfba4d).
  - Snapshot-trigger policy gates (c523cf508).

---

<a id="pr-37"></a>

## PR #37 — Startup smoothing-length depth fix (2026-07-30, merge 6308b70f2)

A single upstream SWIFT fix for startup smoothing-length sanitisation.

### Upstream SWIFT
- When a gas or star particle's `h` is reset at startup (h = 0 or h > upper_h_max), its `depth_h` is now set to the cell depth so the two stay consistent (257b62373).

---

<a id="pr-36"></a>

## PR #36 — Obsidian AGN wind rework, SF-gas FoF galaxies, new FIRE G0 fit (2026-07-30, merge 2ae14470c)

One KIARA commit (1c49330c8) reworks the Obsidian AGN wind speeds and ADAF kicks. It also removes central-BH tracking, restricts FoF galaxy gas to star-forming gas, and replaces the FIRE-derived G0 fit. The rest is an upstream SWIFT sync.

> ⚠️ Restart files from earlier builds are not compatible. The `is_central_bh` field was removed from `bpart`.

### Bug fixes
- **ADAF energy injection could divide by zero** (1c49330c8). Energy is now injected only when `adaf_wt_sum > 0`.
- **Limiter sort-index guard bound was too loose** (1c49330c8). An index is now flagged out of range above `nparts_total/8`, instead of above `nparts_total`.

### New features / options
| parameter | default | what it does |
|---|---|---|
| `ObsidianAGN:adaf_decouple_time_factor` | 0 | Decoupling time for ADAF winds. It was previously hard-wired to 0. |

- **Variable (Simba-like) wind speed now also works for ADAF mode** (1c49330c8). Each mode with a negative wind-speed parameter (`adaf_`, `quasar_`, `slim_disk_wind_speed`) uses the variable M_BH-dependent speed. A positive value is used as a fixed speed. Δlog M_BH is floored at 0. The speed is floored at `minimum_v_kick_km_s` and capped at |`jet_velocity`|. The speed is stored in `bp->v_kick` for all states, including ADAF (it used to be zeroed in ADAF).

### Behaviour changes
- **ADAF kicks are probabilistic** (1c49330c8). A particle is kicked with probability `ObsidianAGN:adaf_kick_factor`, and a kicked particle is not heated. Otherwise it is heated. With a negative `adaf_wind_speed`, the kick speed is min(√(2E/m), |`adaf_wind_speed`|). The old split of energy between kick and heat, and the jet ramp on the ADAF speed, were removed.
- **Cooling shutoff is now `factor × min(t_sound, 10 dt)`, not `factor × min(t_sound, dt)`**, in all AGN modes (1c49330c8).
- **Every particle touched by BH feedback leaves subgrid-ISM mode**, not only kicked ones. Dust is destroyed only in jet particles, no longer in all ADAF-heated gas (1c49330c8).
- **Jet-velocity scaling with BH mass is normalised at 10⁹ M☉** instead of 10⁸ M☉ (`jet_velocity_scaling_with_BH_mass`) (1c49330c8).
- **Quasar/slim-disk decoupling-time assignment removed from the kick loop** (1c49330c8).
- **FoF galaxy gas includes only star-forming gas** (`WITH_FOF_GALAXIES`) (1c49330c8). This applies to the galaxy gas mass, the SFR, the maximum density and the densest-particle position.
- **Central-BH tracking removed** (1c49330c8). `is_central_bh` and the FoF most-massive-BH reduction are gone. `ObsidianAGN:reposition_central_only` is still read but does nothing.
- **New FIRE-derived G0 fit** (`KIARACooling:G0_computation_method = 6`) (1c49330c8). log G0 = log t_ff[Myr] + 3.18·log(log T + Z) − 3.76 + 3.37·exp(−0.173·log v_pot²)·min(1+z, 4). The ad hoc redshift factor and a floor of 100 km²/s² on the potential are new. T is capped at 1000 K and Z floored at 10⁻⁶.

### Other
- Added a `BH_FEEDBACK` debug message for BHs above 10⁸ M☉. The BH log line now includes v_kick (1c49330c8).

### Upstream SWIFT
- Proxy tags exchange (928e48372). DMO runs skip the sort-freeing loop, which is now timed (c049a817b).
- GEAR sinks: bug fixes and clean-up (511a16c90), softening updated on swallowing (692613465), newborn stars use smoothed density/u (ac9b42f05). GEAR added to CI (d7adc1946).
- Averaged output quantities are kept non-negative (e06f6eb31). Compression filters are checked before use (c7faf52c0). Trigger override is no longer recorded on restart (05f01d3d4).
- `cycle.h` updated, now with RISC-V support (9723bd670).
- Examples: idealised-cluster ICs (122e70b55), GEAR/AGORA script clean-up (764567688), Isolatedgalaxy_nfw_mn typo (1d3c892b4), gresho2d plotting (b8ddcc525).
- Contributor Covenant code of conduct adopted (dfc3ae93a). `.gitignore` updates (73902900f, 90030c447).

---

<a id="pr-35"></a>

## PR #35 — Cooling-shutoff follow-ups and BH repositioning by galaxy mass (2026-05-13, merge c88b100dd)

More cooling-shutoff fixes, a galaxy stellar-mass condition on BH repositioning and swallowing, and a workaround for a crash in `part_is_inhibited`.

### Bug fixes
- **Crash from out-of-range sort indices in the limiter loops** (7dfe7285b). New `check_part_sort_index` guard: an out-of-range index is replaced with 0 and a warning is printed. The commit calls this a "kloodgy" workaround, not a root-cause fix.
- **Cooling-shutoff edge cases** (db696078b). `dt` is added to the shutoff timer only if the timer was actually set. The timer is zeroed when a wind particle recouples. Particles with cooling shut off now skip the subgrid-ISM property update, as decoupled particles already did.

### Behaviour changes
- **Slim-disk mode now also shuts off cooling**, using `ObsidianAGN:quasar_cooling_shutoff_factor` (db696078b).
- **`ObsidianAGN:adaf_cooling_shutoff_factor` default changed from -1 to 0** (db696078b).
- **BHs reposition only onto gas in a galaxy of the same or larger stellar mass** (db696078b). A BH is repositioned onto, or swallowed by, another BH only if that BH sits in a galaxy of the same or larger stellar mass (7dfe7285b).
- **Snapshot fields** (7dfe7285b): `SubgridColdISMFraction` renamed `SubgridColdFraction`. The `IonMassFractions` field (non-RT builds) was removed.

### Other
- Debug-message frequency reduced and some warnings toned down (7dfe7285b).

### Upstream SWIFT
- Optional fixed METIS/ParMETIS seed for reproducible debugging (53aa2eebf).
- Corrected the rounding tolerance in the particle-position check for the limiter case (1c72788ec).
- New Balsara–Kim example driven by stars (7b501394a). makeIC scripts updated for the new swiftsimio writer (72659dfda). coolingHaloWithSpin test script fix (57b5d3312).
- Python formatter version bump (06c361ed5).

---

<a id="pr-34"></a>

## PR #34 — AGN cooling shutoff for all modes (2026-05-05, merge 11056dfb3)

Applies the cooling-shutoff timer to all AGN feedback and makes it count down correctly. Also tightens which gas BHs reposition onto and receive feedback from.

### Bug fixes
- **The cooling-shutoff timer only counted down for decoupled particles** (a888206a9). `feedback_recouple_part` now decrements `cooling_shutoff_delay_time` every step for any particle with the timer set, and clamps it at 0. Decoupling and recoupling are handled separately inside the same routine. When a BH sets the timer, `dt` is added because the first decrement happens immediately, so the shutoff lasts at least one step.

### New features / options
- New snapshot field `CoolingShutoffTimes`: time left with radiative cooling and heating switched off (a888206a9).

### Behaviour changes
- **Cooling shutoff by AGN mode** (a888206a9). ADAF-heated gas uses `ObsidianAGN:adaf_cooling_shutoff_factor` and quasar-mode gas uses `ObsidianAGN:quasar_cooling_shutoff_factor`.
- **BH repositioning and feedback skip non-cooling gas** (a888206a9). Repositioning ignores particles that are decoupled or have cooling shut off. It also targets only star-forming gas (SFR > 0). BH feedback is not applied to decoupled or cooling-shutoff particles.
- **`ObsidianAGN:reposition_central_only` no longer has any effect** (a888206a9). The check was removed from the repositioning code, though the parameter is still read.
- **Grackle `H2_self_shielding` changed from 4 to 3** in the KIARA cooling initialisation (a888206a9).
- The firehose chemistry interaction now calls `hydro_diffusive_feedback_reset` on the receiving particle, imposing maximal viscosity (a888206a9).

### Upstream SWIFT
- Added a missing sink condition when deciding whether hydro tasks can be split further (1250eb20c).

---

<a id="pr-33"></a>

## PR #33 — Upstream SWIFT sync (2026-05-04, merge 5e09e5436)

Brings in upstream SWIFT fixes only. No KIARA-specific changes.

### Upstream SWIFT
- Matrix inversion is scaled by the rms so it stays stable when values are small (35228607b).
- `cell_can_use_mesh` is now robust to floating-point errors (0c4e25cf0).
- The power-spectrum fold-factor check runs only when Nfold > 1 (b9beb8f9d).
- Removed duplicated sink logic in `runner_time_integration` (f263b0b05).
- Examples: `pre_supernovae_efficiency` renamed `stellar_winds_efficiency` (0a8ddc8b3); stellar-winds example initial temperature updated (4b9b720d4).
- Documentation link fixes (0e33c4e7c).

---

<a id="pr-32"></a>

## PR #32 — Initial KIARA model import: MAGMA2, KIARA cooling/chemistry/feedback, Obsidian BHs, KIARA-RT (2026-04-26, merge 1b3186832)

First merge of the KIARA galaxy-formation model into master. It rebases the model from the old `romeeld/swiftsim` fork (MAGMA2 hydro, KIARA chemistry/cooling via Grackle/crackle, KIARA stars/star formation/firehose feedback, Obsidian black holes, an experimental KIARA-RT module). It then adds five months of fixes and tuning on top: BH accretion, AGN wind speeds, the G0 (ISRF) models and snapshot I/O.

> ⚠️ Restart files from earlier builds are not compatible. The merge adds new `part`/`xpart`/`spart`/`bpart` fields for the new modules.

### Bug fixes
- **MAGMA2 conduction term differed between symmetric and non-symmetric loops** (8e9c9a209). `du/dt_cond` is now computed the same way in both.
- **Recoupled wind particles were put back onto the timeline a second time** (94d52f705). They are already on it, so the extra re-insertion was removed.
- **u↔T conversion used the wrong mean molecular weight** (e6fa93ca9). `cooling_convert_u_to_temp` and the inverse now compute μ from the tracked species fractions (HI, HII, H⁻, He, H2, e⁻). The old fixed-helium formula is gone, and the inverse was corrected from `T·μ·f` to `T·f/μ`.
- **Snapshot `Temperatures` used the wrong conversion** (eb856a5c1). Output temperatures now come from the physical internal energy and the particle's electron fraction via `cooling_convert_u_to_temp`.
- **`cooling_sputter_dust` skipped hydrogen when looping over elements and mis-accumulated total metallicity** (278384346). The loop now starts at element 0, skips H/He explicitly, and adds each metal to `metal_mass_fraction_total`.
- **BH seeding with FoF galaxies used halo mass instead of stellar mass** (76b4997c2). With `WITH_FOF_GALAXIES`, the seeding threshold is now compared with the group stellar mass.
- **`lum_thresh_always_jet` was applied only at kick time, not to the BH state** (781844de5). A BH above the luminosity threshold now switches on its jet in the state logic. The forced `jet_prob = 1` in the interaction loop was removed.
- **Firehose mixing used `xp->v_full` instead of `p->v`** (781844de5). KIARA chemistry and firehose exchange now use and update `p->v`.
- **Inconsistent "cold gas" definitions in BH accretion** (476cd27d0, c233d33b7, 8e86bdea8). Cold and co-rotating gas masses now use one definition. With the subgrid ISM, only `subgrid_fcold × dense_gas_fraction` of each neighbour counts (it was the H2 fraction). The `cold_disk_mass` field was renamed `corot_gas_mass`.
- **Duplicate `hydro_prepare_gradient` definition in MAGMA2** (d7556e205, c233d33b7). Removed, and the function was restored to its original layout.
- **Configure bug** (9f062f958), plus `subgrid=KIARA` compile fixes, static-function and Obsidian debugging-check fixes, and a build fix for hydro-only configurations (e5dabe62e, 008761e88, 9b1999be9, e40834248, 9726edee4, 2ffd4e01d, f4637aa60).

### New features / options
- **New subgrid physics modules**: MAGMA2 SPH (`--with-hydro=magma2`), KIARA chemistry (firehose wind–ISM mixing), KIARA cooling (Grackle/crackle with subgrid ISM, dust, H2 and G0), KIARA stars and star formation, KIARA feedback (decoupled winds with recoupling), Obsidian black holes, and KIARA-RT (f0e6d35c2, ac97b179b, e43a18d7c, 632fe0287, d03df6f9e, 0e8f3027d).
- **New configure options**: `--with-subgrid=SIMBA|KIARA|KIARART`, `--with-subgrid-black-holes=Obsidian` presets (also `ObsidianNoBH`), `--enable-fof-galaxies` (FoF galaxy finding, `WITH_FOF_GALAXIES`), and `--enable-{simba,kiara,firehose,obsidian,magma2}-debugging-checks`. Grackle is now linked with an rpath.
- **New run-time flags**: `--kiara` and `--kiarart` switch on hydro, limiter, sync, self-gravity, stars, star formation, cooling, feedback, black holes and FoF. The error that stopped RT and cooling from running together was removed.
- **New snapshot outputs**: `Temperatures`, G0 values (c83d9e0ec, e6fa93ca9), and unit-tagged subgrid temperature, dust temperature and density fields. `ElectronNumberDensities` is now given relative to n_H (d7556e205).
- **KIARA-RT** (7504a32d4, eb856a5c1) can use native KIARA/Grackle cooling. It works without subcycling and fails with a sorting error with subcycling. Radiation energy densities and fluxes are now corrected for changes in smoothing length.
- **Obsidian galaxy-centre tracking** (8e86bdea8, d68eb1005). FoF galaxies track their most massive BH (`is_central_bh`). Torque accretion is switched off when the BH is far from the centre (9f062f958).
- **Variable quasar/slim-disk wind speed** (95bca10b0, 278384346). A negative `quasar_wind_speed` / `slim_disk_wind_speed` gives v = |v_q| + (|v_sd|/3)·Δlog M_BH. A negative `slim_disk_wind_speed` makes the Δlog M_BH term quadratic. The speed is capped at the jet velocity.
- **FIRE-based G0 fit** (c83d9e0ec). `G0_computation_method = 6` uses a symbolic-regression fit in the free-fall time and the gravitational potential (marked preliminary).

Options added during the development cycle on top of the imported model:

| parameter | default | what it does |
|---|---|---|
| `SPH:target_neighbours` | 0 (unset) | Sets the neighbour number directly. Either this or `SPH:resolution_eta` must be given, and `SPH:resolution_eta` is now optional. MAGMA2 requires the compiled `const_kernel_target_neighbours` to be within 5% of this value. |
| `ObsidianAGN:bondi_BH_mass_cap_Msun` | FLT_MAX | Caps the BH mass used in the Bondi rate. |
| `ObsidianAGN:cold_gas_temperature_cut_K` | = `environment_temperature_cut_K` | Temperature below which gas counts as cold for accretion. |
| `ObsidianAGN:quasar_cooling_shutoff_factor` | = `adaf_cooling_shutoff_factor` | Cooling-shutoff multiplier for gas heated in quasar mode. |
| `ObsidianAGN:reposition_central_only` | 0 | Repositions only the central (most massive) BH. |
| `KIARACooling:G0_multiplier` | 1.0 | Scales the computed G0. |
| `KIARACooling:G0_computation_method` | 3 | Now read as an int. New methods `-3` and `6` (FIRE fit). |
| `KIARACooling:do_cooling_in_rt` | 0 | Runs cooling inside RT thermochemistry (KIARA-RT). |

### Behaviour changes
- **BH accretion uses the dense gas fraction, not the H2 fraction** (8e86bdea8). The commit calls this a major change.
- **Torque-accretion suppression refactored** (c233d33b7). A single `f_suppress` factor is now applied. A negative `ObsidianAGN:torque_accretion_norm` adds an exponential cutoff above M_BH = 10⁹ M☉.
- **Gas heated by AGN stops cooling** for `adaf_cooling_shutoff_factor × min(t_sound, dt)` so it does not drop straight back into ISM mode (8e86bdea8).
- **`ObsidianAGN:minimum_v_kick_km_s` default raised from 10 to 100** (278384346).
- **`KIARAFeedback:FIRE_eta_lower_slope_EOR` removed** (c78aa5f22).
- **Timestep-limiter "not synced" check disabled** (commented out) in `runner_time_integration.c`. Wind recoupling now runs in the time-integration loop (`feedback_recouple_part`).

### Other
- Integration and clean-up by Matthieu Schaller: snapshot I/O names follow SWIFT conventions (`CoolingTime`→`CoolingTimes`, 7532928e1, 09da22988). Expensive gradient maths moved from `end_density` to `prepare_gradient` (5affc77b6). Symmetry/MPI-rules unit tests now cover KIARA chemistry (bb617ca3c, 86bc2ec4b, aa004612e). Also `hydro_debug` decoupling info (af339d1b0), `.gitignore` for KIARA side tables (9416fa3b5), code formatting (5e9e9d916, 2f1e4defc), and stray files removed (30f586939).
- Stars ghost: bisection fallback added for when Newton–Raphson over-iterates on h (d7556e205). MAGMA2 `div_v` added for unit tests (c577ce601). Empty star-formation logger routines added (d68eb1005).
- Diagnostics: extended BH accretion log line (galaxy M*, SFR, h, f_suppress, t_dyn), and periodic G0 debug messages.

