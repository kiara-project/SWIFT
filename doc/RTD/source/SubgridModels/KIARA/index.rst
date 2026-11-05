.. KIARA size-resolved dust
   Balu Sreedhar, 2026

.. _KIARA_dust:

KIARA dust
==========

The KIARA dust model evolves the dust mass of each element on the gas
particles with the built-in crackle chemistry. Optionally, each gas particle
also carries the size distribution of its carbonaceous and silicate grains on
a common grid of grain radii. The size bins describe how each material's dust
mass is distributed in radius; the dust mass itself, and its elemental
composition, remain those evolved by crackle.

*The size-resolved dust is under development: the grid and the stellar
source spectra below are in place, their evolution in the gas is being
added.*

Grain-size grid
---------------

Choose the number of radii per material at build time with
``--with-kiara-dust-size-bins=N``. The default, 0, leaves the size-resolved
dust out of the build (no memory is used for it); otherwise ``N`` must be at
least 2. The number sets the particle layout, so changing it requires
rebuilding and restart files cannot be converted.

Set ``KIARACooling:dust_grain_size_min_micron`` and
``KIARACooling:dust_grain_size_max_micron`` at run time to construct
logarithmically spaced representative radii including both endpoints. The
defaults are 0.005 and 0.1 micron. For the wider range of Li et al. (2021),
set these to 0.0001 and 1.

Alternatively, ``KIARACooling:dust_grain_sizes_micron`` supplies exactly
``N`` finite, positive, strictly increasing radii and overrides the endpoint
keys. An explicit grid may be non-uniform. The radii are physical grain
radii, not comoving lengths, and are shared by the cooling and the stellar
dust injection.

Surviving dust below the minimum represented radius is deposited in the
smallest bin; dust above the maximum is deposited in the largest bin. This
retains mass at the boundaries but cannot retain the unresolved grain number
or surface area. Physically destroyed dust is returned to gas-phase elements.
Six bins across a wider range remain a coarse approximation; a configurable
geometry does not imply convergence with respect to resolution.

``KIARACooling:dust_grainsize`` is the separate representative radius used by
crackle and does not define the size grid.

Dust sources
------------

``KIARACooling:dust_initial_source`` selects the source distribution used for
dust without an existing size distribution, for empty material rows, and for
the artificial self-enrichment seed. It defaults to ``SNII`` and also accepts
``AGB``. This is an explicit assumption about dust without a recorded source
history, not an inference of that history from its composition. The source
parameters and integration are shared with the stellar injection.

Stellar injection integrates the mass spectrum corresponding to Li et al.
(2021), equation 5, over geometric bin cells. The decaying lognormal form is
used (the printed equation omits the minus sign in the exponent). Both tails
are added to the endpoint bins rather than discarded or sampled randomly. The
number spectrum is proportional to
``a**(-p) * exp(-ln(a/a0)**2/(2*sigma**2))`` per logarithmic radius interval;
mass weighting contributes another factor of ``a**3``.

The source parameters are ``KIARAFeedback:dust_SNII_a0_micron`` (0.1),
``dust_SNII_sigma`` (0.6) and ``dust_SNII_p`` (0), with the corresponding
``dust_AGB_a0_micron`` (0.1), ``dust_AGB_sigma`` (0.47) and ``dust_AGB_p``
(4). Both materials use the same source spectrum for a given channel; their
masses and subsequent size distributions remain independent. The 100 Myr
stellar-age split between the SNII and AGB channels is an approximation.

Validation
----------

``make check`` runs the stand-alone checks of the grid, the source spectra
and the radius transport and grain collisions with 6 bins. Run
``python3 tests/testKIARADustGrid.py`` to compile and run them with 2, 3, 6,
12 and 41 bins under the address and undefined-behaviour sanitizers. These
do not link crackle or GSL and do not replace a full KIARA build or a
convergence study.
