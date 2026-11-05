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

Grain-size grid
---------------

Choose the number of radii per material at build time with
``--with-kiara-dust-size-bins=N``. The default is 10; 0 leaves the
size-resolved dust out of the build (no memory is used for it, and crackle
uses the single grain size ``KIARACooling:dust_grainsize``); otherwise ``N``
must be at least 2. Each bin adds 16 bytes to a gas particle (10 bins: 1184
instead of 1024 bytes) and nothing to a star particle. The number sets the
particle layout, so changing it requires rebuilding and restart files cannot
be converted.

Set ``KIARACooling:dust_grain_size_min_micron`` and
``KIARACooling:dust_grain_size_max_micron`` at run time to construct
logarithmically spaced representative radii including both endpoints. The
defaults are 0.001 and 1 micron (Li et al. 2021 use 0.0001 to 1 micron).

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

Evolution in crackle
--------------------

With size bins, crackle evolves the dust from its size distribution (there is
no option to keep the bins as bookkeeping only). In Li et al. (2021),
accretion growth and thermal sputtering change the radius of every grain at
the same rate, which crackle's rates give for 0.1 micron grains:

* growth: each element's accretion time scales with the area-weighted radius
  ``a_eff = sum(m) / sum(m / a)`` of its material (carbon in carbonaceous
  grains; O, Mg, Si, S, Ca and Fe in silicates and oxides), and the surviving
  grains then grow by the common ``da`` that adds the accreted mass;
* sputtering: every grain is eroded by the same
  ``da = 0.1 micron dt / t_sput(0.1 micron)``; grains smaller than ``da`` are
  destroyed;
* SN shocks: every grain loses the fraction of its mass that the
  (size-independent) shock destruction removes.

The mass destroyed in each material is taken from its elements in proportion
to their dust masses, and after crackle's per-element limits each material's
bins are rescaled to its dust mass, so that the bins always agree with the
element composition. Thermal sputtering outside crackle (hot and Firehose
gas) works in the same way.

The surface processes scale with the grain area per dust mass relative to 0.1
micron grains, ``0.1 micron / a_eff`` over both materials: H2 formation on
grains, the gas-grain heat exchange (and so the dust temperature),
photoelectric heating and recombination onto grains (and the Sternberg et
al. 2014 limit on the molecular fraction through its H2 formation rate).
Absorption of the interstellar radiation field (dust heating, dust shielding
of the Lyman-Werner band) is per dust mass, as for grains much smaller than
the wavelength. ``KIARACooling:dust_grainsize`` is not used with size bins.

Grains cannot grow beyond the largest radius of the grid, so its upper end
should allow for growth.

Grain collisions
----------------

Grain-grain collisions move dust between the bins without changing its mass
(Li et al. 2021, section 2.5.3). The grains are coupled to the turbulence of
the gas, with a velocity dispersion (Li et al., equation 15)

.. math::
   \sigma_{\rm gr} = 0.06\,{\rm km\,s^{-1}} \left(\frac{v_J}{c_s}\right)^2
   \frac{a}{0.1\,\mu{\rm m}} \left(\frac{n}{1\,{\rm cm^{-3}}}\right)^{-1/2}
   \frac{\rho_{\rm gr}}{2.4\,{\rm g\,cm^{-3}}},

where :math:`v_J = 0.7\,(L_J/1\,{\rm pc})^{1/2}` km/s is the turbulent
velocity at the Jeans length :math:`L_J`. Every step, each bin gets a velocity
whose components are drawn from a Gaussian of dispersion
:math:`\sigma_{\rm gr}/\sqrt{3}` (reproducibly, from the particle ID and the
time). Pairs of grains slower than the sticking threshold of Hirashita & Yan
(2009, equation 8) coagulate; pairs faster than 2 km/s shatter, into
fragments with a :math:`a^{-3.3}` spectrum following Hirashita & Yan (2009),
section 2.3. The material constants are those of Hirashita & Yan's table 2
(quartz for silicates, graphite for carbonaceous grains), with the common
grain density ``KIARACooling:dust_grain_density`` (2.4 g/cm^3). Grains only
collide with grains of the same material.

The cascade of Larson's law up to the Jeans length only makes sense in the
ISM: in diffuse gas the Jeans length is many kpc and the formula gives grains
in the circumgalactic and intergalactic medium hundreds of km/s, which grinds
them below the smallest bin. So the collisions only happen above the density
threshold of the subgrid ISM (``KIARACooling:subgrid_threshold_n_H_cgs``, with
more than 100 times the mean baryon density), and the grain velocities are at
most ``KIARACooling:dust_collision_velocity_cap`` (1) times the velocity
dispersion of the gas across the kernel, :math:`H|S|` from the trace-free
velocity shear :math:`S` of the Smagorinsky diffusion (no cap if diffusion is
off; 0 also switches the cap off). In the subgrid ISM the collisions happen
in the cold phase, at the density and temperature that crackle uses, where
the grains are slow and mostly coagulate; shattering happens in the warm and
hot ISM.

``KIARACooling:dust_coagulation`` and ``KIARACooling:dust_shattering`` (both 1
by default) switch the two processes.

Transport
---------

Wherever dust moves, its sizes move with it:

* stellar dust is added to the gas with the spectrum of its channel (SNII for
  stars younger than 100 Myr, AGB otherwise);
* diffusion and Firehose mixing exchange the dust of each bin with the same
  pairwise weights as the dust of each element;
* self-enrichment seeds dust with the ``dust_initial_source`` spectrum, and
  jets destroy all of a particle's dust.

In each case the moved dust sets the shape of the distribution, while the
dust masses of the elements set the amount of each material, so the bins and
the element composition always agree.

Validation
----------

``make check`` runs the stand-alone checks with 6 bins: the grid, the source
spectra, the radius transport (growth, sputtering), the bins kept consistent
with the element dust, the SN shock shrink, and the grain collisions (Li et
al. Eq. 15 and Hirashita & Yan Eq. 8 against their published forms, single
pair rates, conservation of each material's mass, coagulation moving mass to
larger grains and shattering to smaller ones). Run
``python3 tests/testKIARADustGrid.py`` to compile and run them with 2, 3, 6,
12 and 41 bins under the address and undefined-behaviour sanitizers. These
do not link crackle or GSL and do not replace a full KIARA build or a
convergence study.
