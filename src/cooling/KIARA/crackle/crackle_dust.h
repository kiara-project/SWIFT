
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DUST_FRAC_CAP 0.9 // Total dust growth is limited to this fraction of gas-phase metals
#define DUST_CHANGE_HI 100.0  // factor by which dust density can go up in a single step
#define DUST_CHANGE_LO 0.5  // factor by which dust density can go down in a single step
#define T_SPUTTERING 2.e4  // Temeprature below which sputtering is negligible
#define DUST_GROWTH_TAU_GYR 1.0  // dust growth time-scale (Gyr) at density dust_growth_densref

/* Apply the changes drho of the dust metal densities (by crackle slot), limited per step and capped at a fraction of
 * the metals, move the metals between gas and dust, and apply the oxygen stoichiometric cap. */
static inline void apply_dust_changes(grackle_part_data *gp, chemistry_data *chemistry, double drho[NUM_METAL_SPECIES_GRACKLE])
{
        int k;
	/* Now evolve the dust metal density */
	gp->dust_density = 0.;
	double old_dustDensity, old_metalDensity;
        for (k=1; k<NUM_METAL_SPECIES_GRACKLE; k++){
	    old_dustDensity = gp->dust_metalDensity[k];
	    old_metalDensity = gp->dust_metalDensity[k] + gp->gas_metalDensity[k];
	    /* Add to dust metallicity */
	    gp->dust_metalDensity[k] += drho[k];
	    /*if (gp->dust_metalDensity[k] < tiny) {
		    drho[k] = - gp->dust_metalDensity[k];
		    gp->dust_metalDensity[k] = tiny;
	    }*/
	    /* Limit change to dust within a single step, with cap on overall fraction of metals in dust */
	    if (gp->dust_metalDensity[k] > DUST_CHANGE_HI*old_dustDensity || gp->dust_metalDensity[k] > DUST_FRAC_CAP*old_metalDensity) {
	    	gp->dust_metalDensity[k] = min(DUST_CHANGE_HI*old_dustDensity, DUST_FRAC_CAP*old_metalDensity);
		drho[k] = gp->dust_metalDensity[k] - old_dustDensity;
	    }
	    if (gp->dust_metalDensity[k] < old_dustDensity*DUST_CHANGE_LO) {
	    	gp->dust_metalDensity[k] = old_dustDensity*DUST_CHANGE_LO;
		drho[k] = gp->dust_metalDensity[k] - old_dustDensity;
	    }
	    /* Metals going into dust are taken out of gas */
	    gp->gas_metalDensity[k] -= drho[k];
	    gp->metal_density -= drho[k];
	    /* Add up new dust densities */
	    gp->dust_density += gp->dust_metalDensity[k];
        }

	/* Grains hold oxygen only bound to refractory elements (silicates and
	 * oxides), so limit the dust O to the oxides of the dust Mg, Si, Ca and Fe
	 * (MgO, SiO2, CaO, FeO) and return any excess to the gas. Indices are those
	 * of gas/dust_metalDensity: 3 = O, 5 = Mg, 6 = Si, 8 = Ca, 9 = Fe. */
	if (chemistry->dust_oxygen_stoichiometric_cap) {
	    const double O_dust_max = 15.999 * (gp->dust_metalDensity[5] / 24.305 +
	                                        2. * gp->dust_metalDensity[6] / 28.085 +
	                                        gp->dust_metalDensity[8] / 40.078 +
	                                        gp->dust_metalDensity[9] / 55.845);
	    const double O_excess = gp->dust_metalDensity[3] - O_dust_max;
	    if (O_excess > 0.) {
	        gp->dust_metalDensity[3] -= O_excess;
	        gp->gas_metalDensity[3] += O_excess;
	        gp->metal_density += O_excess;
	        gp->dust_density -= O_excess;
	    }
	}
}

#if KIARA_DUST_N_BINS > 0
#include "chemistry/KIARA/dust_evolution.h"

/* Dust material of a crackle metal slot: carbon (slot 1) forms carbonaceous grains, the other dust elements (O, Mg, Si,
 * S, Ca, Fe) silicates and oxides */
static inline int crackle_dust_material(int k) { return k == 1 ? kiara_dust_carbonaceous : kiara_dust_silicate; }

/* Dust evolution with size-resolved grains (KIARA_DUST_N_BINS > 0): as evolve_dust(), with every process acting on the
 * grain-size distribution of each material (gp->dust_bins, fractions of the dust mass in bins of radius
 * gp->dust_grain_sizes, in micron). In Li et al. (2021) accretion growth and thermal sputtering change the radius of
 * every grain at the same rate, da/dt, which crackle's rates give for a_ref = 0.1 micron grains:
 *  - growth: each element's accretion time scales with the area-weighted radius a_eff = sum(m) / sum(m / a) of its
 *    material (t_accr ~ a), and the surviving grains then grow by the common da that adds the accreted mass;
 *  - sputtering: every grain is eroded by the same da = a_ref dt / t_sput(a_ref); grains smaller than da are destroyed;
 *  - SN shocks: every grain loses the fraction of its mass that crackle's (size-independent) destruction removes.
 * The mass destroyed in each material comes from its bins and is taken from its elements in proportion to their dust
 * masses. After the per-element limits of apply_dust_changes(), each material's bins are rescaled to its dust mass, so
 * the bins and the element masses always agree. */
static inline void evolve_dust_binned(grackle_part_data *gp, chemistry_data *chemistry, code_units *units, int ism_flag, double dtit)
{
	const double sec_per_year = 3.1536e7;
	const double a_ref = 0.1;  // micron: crackle's growth and sputtering rates are those of grains of this radius
	const double *radius = gp->dust_grain_sizes;
	double drho[NUM_METAL_SPECIES_GRACKLE], drho_grow[NUM_METAL_SPECIES_GRACKLE];
	int k, mat, i;

	if (gp->dust2gas <= tiny || gp->dust_density <= tiny) return;
	const double dens_cgs = gp->density * units->density_units;
	for (k=0; k<NUM_METAL_SPECIES_GRACKLE; k++) drho[k] = drho_grow[k] = 0.;

	/* Dust (density) in each material and in its bins; the bins hold fractions of the total dust mass */
	double m_mat[KIARA_DUST_N_MATERIALS], mass[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS];
	double after[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS], out[KIARA_DUST_N_BINS];
	for (mat=0; mat<KIARA_DUST_N_MATERIALS; mat++) m_mat[mat] = 0.;
	for (k=1; k<NUM_METAL_SPECIES_GRACKLE; k++) m_mat[crackle_dust_material(k)] += gp->dust_metalDensity[k];
	for (mat=0; mat<KIARA_DUST_N_MATERIALS; mat++) {
	    double row = 0.;
	    for (i=0; i<KIARA_DUST_N_BINS; i++) { row += gp->dust_bins[mat][i]; mass[mat][i] = 0.; }
	    if (row > 0.) {
		for (i=0; i<KIARA_DUST_N_BINS; i++) mass[mat][i] = gp->dust_bins[mat][i] * m_mat[mat] / row;
	    } else if (m_mat[mat] > 0.) {
		kiara_dust_grid_deposit(radius, a_ref, m_mat[mat], mass[mat]);  // no shape known: grains of radius a_ref
	    }
	}

	/* Area-weighted radius of each material */
	double a_eff[KIARA_DUST_N_MATERIALS];
	for (mat=0; mat<KIARA_DUST_N_MATERIALS; mat++) {
	    double moment = 0.;
	    for (i=0; i<KIARA_DUST_N_BINS; i++) moment += mass[mat][i] / radius[i];
	    a_eff[mat] = (moment > 0.) ? m_mat[mat] / moment : a_ref;
	}

	/* Dust growth: metals accreted from gas to dust; only allowed when in ISM. As in evolve_dust(), for grains of
	 * radius a_eff */
	if (ism_flag) {
	    const double tau_ref = DUST_GROWTH_TAU_GYR * 1.e9 * sec_per_year / units->time_units;
	    const double tau_accr0 = tau_ref * (chemistry->dust_growth_densref/dens_cgs) * sqrt(gp->tdust/gp->tgas);
	    for (k=1; k<NUM_METAL_SPECIES_GRACKLE; k++){
		const double rhomet = gp->gas_metalDensity[k] + gp->dust_metalDensity[k];
		const double tau_accr = tau_accr0 * chemistry->SolarAbundances[k] * gp->density / (gp->gas_metalDensity[k]+tiny)
		                        * (a_eff[crackle_dust_material(k)] / a_ref);
		drho_grow[k] = min((gp->gas_metalDensity[k] / rhomet) * (gp->dust_metalDensity[k]/tau_accr) * dtit, DUST_FRAC_CAP * gp->gas_metalDensity[k]);
	    }
	}

	/* Dust destruction, acting on the bins */
	for (mat=0; mat<KIARA_DUST_N_MATERIALS; mat++)
	    for (i=0; i<KIARA_DUST_N_BINS; i++) after[mat][i] = mass[mat][i];
	/* SNe shocks: the fraction of the dust that evolve_dust() destroys, as the same fractional loss of every grain */
	if (gp->SNe_density > tiny) {
	    const double mass_unit = units->density_units * units->length_units * units->length_units * units->length_units;
	    const double rhogas = gp->density - gp->dust_density;
	    const double Ms100 = 6800.0*chemistry->sne_coeff*(100.0/chemistry->sne_shockspeed)*(100.0/chemistry->sne_shockspeed) * SolarMass / mass_unit;
	    const double f_shocked = min((Ms100 * gp->SNe_density) / rhogas, 1.);
	    const double f_lost = min(f_shocked * chemistry->dust_destruction_eff, 1.);
	    const double scale = cbrt(1. - f_lost);
	    for (mat=0; mat<KIARA_DUST_N_MATERIALS; mat++) {
		kiara_dust_radius_remap(radius, after[mat], 0., scale, 1, out);
		for (i=0; i<KIARA_DUST_N_BINS; i++) after[mat][i] = out[i];
	    }
	}
	/* Sputtering (Tsai & Mathews 1995): every grain is eroded by the same da; negligible below T_SPUTTERING in the ISM */
	if (gp->tgas > T_SPUTTERING || !ism_flag) {
	    const double tau_sput_ref = 1.7e8 * sec_per_year / units->time_units * (1.e-27 / dens_cgs)
	                                * (pow(2.e6 / gp->tgas, 2.5) + 1.0);  // for grains of radius a_ref
	    const double da = a_ref * dtit / tau_sput_ref;
	    for (mat=0; mat<KIARA_DUST_N_MATERIALS; mat++) {
		kiara_dust_radius_remap(radius, after[mat], -da, 1., 1, out);
		for (i=0; i<KIARA_DUST_N_BINS; i++) after[mat][i] = out[i];
	    }
	}
	/* The destroyed mass of each material is taken from its elements in proportion to their dust masses */
	double dm_dest[KIARA_DUST_N_MATERIALS];
	for (mat=0; mat<KIARA_DUST_N_MATERIALS; mat++) {
	    double left = 0.;
	    for (i=0; i<KIARA_DUST_N_BINS; i++) left += after[mat][i];
	    dm_dest[mat] = fmax(m_mat[mat] - left, 0.);
	}
	for (k=1; k<NUM_METAL_SPECIES_GRACKLE; k++) {
	    mat = crackle_dust_material(k);
	    drho[k] = drho_grow[k];
	    if (m_mat[mat] > 0.) drho[k] -= gp->dust_metalDensity[k] / m_mat[mat] * dm_dest[mat];
	}

	/* The surviving grains grow by the common da that adds the accreted mass */
	for (mat=0; mat<KIARA_DUST_N_MATERIALS; mat++) {
	    double dm_grow = 0.;
	    for (k=1; k<NUM_METAL_SPECIES_GRACKLE; k++) if (crackle_dust_material(k) == mat) dm_grow += drho_grow[k];
	    const double da = kiara_dust_growth_shift(radius, after[mat], dm_grow);
	    if (da > 0.) {
		kiara_dust_radius_remap(radius, after[mat], da, 1., 1, out);
		for (i=0; i<KIARA_DUST_N_BINS; i++) after[mat][i] = out[i];
	    }
	}

	apply_dust_changes(gp, chemistry, drho);

	/* The bins: the evolved shapes, rescaled to each material's dust after the limits of apply_dust_changes() */
	double m_new[KIARA_DUST_N_MATERIALS];
	for (mat=0; mat<KIARA_DUST_N_MATERIALS; mat++) m_new[mat] = 0.;
	for (k=1; k<NUM_METAL_SPECIES_GRACKLE; k++) m_new[crackle_dust_material(k)] += gp->dust_metalDensity[k];
	for (mat=0; mat<KIARA_DUST_N_MATERIALS; mat++) {
	    double left = 0.;
	    for (i=0; i<KIARA_DUST_N_BINS; i++) { left += after[mat][i]; gp->dust_bins[mat][i] = 0.; }
	    if (gp->dust_density <= tiny || m_new[mat] <= 0.) continue;
	    if (left > 0.) {
		for (i=0; i<KIARA_DUST_N_BINS; i++) gp->dust_bins[mat][i] = after[mat][i] * (m_new[mat] / left) / gp->dust_density;
	    } else {
		kiara_dust_grid_deposit(radius, a_ref, m_new[mat] / gp->dust_density, gp->dust_bins[mat]);
	    }
	}
	gp->dust_area_factor = kiara_dust_area_factor(radius, (const double (*)[KIARA_DUST_N_BINS])gp->dust_bins, a_ref);
}
#endif

static inline void evolve_dust(grackle_part_data *gp, chemistry_data *chemistry, code_units *units, int ism_flag, double dtit)
{
#if KIARA_DUST_N_BINS > 0
	evolve_dust_binned(gp, chemistry, units, ism_flag, dtit);
	return;
#endif
        // dust growth via accretion and destruction via shock and thermal sputtering
        int k;
        double drho[NUM_METAL_SPECIES_GRACKLE];  // change in dust metal density
	double rhomet[NUM_METAL_SPECIES_GRACKLE];  // total metallicity including both gas and dust
	double drhos;  // dust mass loss via shock and sputtering
        double tau_ref=1.e20, tau_accr=1.e20, tau_accr0=1.e20, f_shocked=0., tau_sput=1.e20;
	double sec_per_year = 3.1536e7;

	if (gp->dust2gas <= tiny || gp->dust_density <= tiny) return;

        //double dens_cgs = gp->density * units->density_units / (units->a_value * units->a_value * units->a_value);  // to physical cgs
        double dens_cgs = gp->density * units->density_units;

        for (k=0; k<NUM_METAL_SPECIES_GRACKLE; k++) drho[k] = 0.;

	// Dust growth: metals accreted from gas to dust; only allowed when in ISM
	if (ism_flag) {
            tau_ref = DUST_GROWTH_TAU_GYR * 1.e9 * sec_per_year / units->time_units;
            tau_accr0 = tau_ref * (chemistry->dust_growth_densref/dens_cgs) * sqrt(gp->tdust/gp->tgas); 
            for (k=1; k<NUM_METAL_SPECIES_GRACKLE; k++){
		//gp->dust_metalDensity[k] = fmax(gp->dust_metalDensity[k], tiny);
                rhomet[k] = gp->gas_metalDensity[k] + gp->dust_metalDensity[k];
                tau_accr = tau_accr0 * chemistry->SolarAbundances[k] * gp->density / (gp->gas_metalDensity[k]+tiny);
		//if (rhomet[k] < tiny || tau_accr < tiny) continue;
                drho[k] = min((gp->gas_metalDensity[k] / rhomet[k]) * (gp->dust_metalDensity[k]/tau_accr) * dtit, DUST_FRAC_CAP * gp->gas_metalDensity[k]); // growth in a single step capped at fraction of the available gas metals
            }
	}

	/* Dust destruction */
	if (gp->dust_density > tiny) {  // can't destroy if there is no dust
	    /* SNe shocks */
            if (gp->SNe_density <= tiny) drhos = 0.0;  // Total change in density from all destruction
            else{
	        const double mass_unit = units->density_units * units->length_units * units->length_units * units->length_units;
                const double rhogas = gp->density - gp->dust_density;
		// gas mass shocked above sne_shockspeed per SN in the Sedov-Taylor phase, ~E/v_s^2 and so independent of the ambient density (McKee 1989; as in Grackle and Simba), in code units
		const double Ms100 = 6800.0*chemistry->sne_coeff*(100.0/chemistry->sne_shockspeed)*(100.0/chemistry->sne_shockspeed) * SolarMass / mass_unit;
                f_shocked = min((Ms100 * gp->SNe_density) / rhogas, 1.);  // fraction of mass shock-heated
                drhos = gp->dust_density * f_shocked * chemistry->dust_destruction_eff;  // some fraction of dust is destroyed in shocked gas
	//if (gp->verbose) printf("dust: %g %g %g %g %g %g\n",gp->dust_density, Ms100 * mass_unit / SolarMass, dens_cgs, f_shocked, drhos, (gp->metal_density/(gp->metal_density+gp->dust_density)) * (gp->dust_density/tau_accr) * dtit);
	    }
	    /* sputtering */
	    if (gp->tgas > T_SPUTTERING || !ism_flag) {  // negligibly small below this temperature
	        tau_sput = 1.7e8 * sec_per_year / units->time_units
                           * (chemistry->dust_grainsize / 0.1) * (1.e-27 / dens_cgs)
                           * (pow(2.e6 / gp->tgas, 2.5) + 1.0); // sputtering timescale, Tsai & Mathews (1995)
                //if (tau_sput > 0.) drhos += gp->dust_density / tau_sput * 3.0 * dtit;
                drhos += gp->dust_density * (1. - exp(-3.* fmin(dtit / tau_sput, 5.)));
	    }
	    /* Destruction is applied to all metals, in proportion to their metallicity */
	    for (k = 1; k < NUM_METAL_SPECIES_GRACKLE; k++) {
                drho[k] -= gp->dust_metalDensity[k] / gp->dust_density * drhos;
	    }
	}

	apply_dust_changes(gp, chemistry, drho);

	//if (gp->verbose) printf("dust: dtit=%g tau_accr=%g tau_sput=%g snerho=%g drhos=%g rhog=%g td=%g tg=%g dust=%g\n",dtit,tau_accr0,tau_sput/3.,gp->SNe_density,drhos,gp->density,gp->tdust,gp->tgas,gp->dust_density);
	//assert(gp->dust_density == gp->dust_density);
	//assert(gp->metal_density >= 0.);
	//assert(gp->dust_metalDensity[0] == 0.); // He should not participate in dust (nor N, Ne)
	//assert(gp->dust_density < gp->density);
	//assert(gp->verbose == 0);

	return;
}

static inline double dust_thermal_balance(double tdust, double tgas, double tcmb, double tcmb4, double gamma_isrf, double gasgr, double nh)
{
        const double kgr1 = 0.0428266; /* dust opacity normalized to local DGR =4.e-4/0.00934 */
        const double kgr200 = 16.f; /* dust-rad coupling at sublimation temp */
        const double tsubl = 1500.f;  /* dust sublimation temperature */
        const double sigma_sb = 5.670373e-5;  /* Stefan-Boltzmann const, in erg/s/cm^2/K^4  */
        double kgr, sol, tdust4;

        /* Initialize: If input tdust<=0, return initial guess for tdust */
        if (tdust <= 0.f) {
            tdust = fmax(tcmb, pow(0.25 * gamma_isrf / sigma_sb / kgr1, 0.17f));
            return tdust;
        }

        /* Otherwise compute dust thermal balance, return heating - cooling rate */
        /* Compute dust-radiation coupling constant */
        if (tdust < 200.f) kgr = kgr1 * tdust * tdust;
        else if (tdust > tsubl) kgr = kgr200;
        else {
            const double x = tdust/tsubl, x4 = x*x*x*x;
            kgr = fmax(tiny, kgr200 * x4*x4*x4);  // (tdust/tsubl)^12
        }
        /* Compute dust heating-cooling rate */
        tdust4 = tdust * tdust * tdust * tdust;
        sol = gamma_isrf + 4.f * sigma_sb * kgr * (tcmb4 - tdust4) + gasgr * nh * (tgas-tdust);
	return sol;
}

static inline void calculate_tdust_bisect(grackle_part_data *gp, double gasgr, double gamma_isrf, double tcmb)
		//double tgas, double nh, double gasgr, double gamma_isrf, double tcmb, double td)
{
        /* Solve for Tdust from ISRF heating, CMB heating, and gas-grain cooling */
        const float tol = 1.e-2;
        const float range = 2.0;   // max change in Tdust (either X or /) allowed in this function call
        const float tcmb4 = tcmb * tcmb * tcmb * tcmb;
        float tdold = -1.e20;
	float td = gp->tdust;
        int iter = 0, maxiter = 100;
        double sol;

	/* Find gas temperature of birth cloud where we assume dust lives, assuming pressure equil with ISM gas */
	const float nbirthcloud = 1.e5;  // This should be above the threshold density in a lognormal SF model.
	float tbirthcloud = gp->tgas * gp->nH / nbirthcloud;
	if (tbirthcloud > gp->tgas) tbirthcloud = gp->tgas;
	if (tbirthcloud < tcmb) tbirthcloud = tcmb;
        /* If first time, initalize Tdust to some reasonable value */
        /* (the bisection range below is set from gp->tdust, so store the guess there) */
        if (gp->tdust <= 0.f) {
            td = dust_thermal_balance(td, tbirthcloud, tcmb, tcmb4, gamma_isrf, gasgr, nbirthcloud);
            gp->tdust = td;
        }
	/* If Tgas has dropped below Tcmb (though it shouldn't); set Tdust to Tgas */
	if (gp->tgas < tcmb) {
	    gp->tdust = gp->tgas;
	    return;
	}
	/* Limit Tdust within allowed range */
	if (gp->tdust > gp->tgas) gp->tdust = gp->tgas;
	if (gp->tdust < tcmb) gp->tdust = tcmb;

        /* Set bisection limits */
	float tdhi = min(gp->tgas, range * gp->tdust);
	float tdlo = max(tcmb, gp->tdust / range);
	/* In rare cases Tdust change is so large that it is outside of allowed range; just return Tbirthcloud */
	if (tdlo >= tdhi) {
	    gp->tdust = tbirthcloud;
	    return;
	}
        /* Solve for Tdust via bisection */
        while (fabs(td-tdold) > tol * td) {
            tdold = td;
	    td = 0.5 * (tdlo + tdhi);  // new guess
            sol = dust_thermal_balance(td, tbirthcloud, tcmb, tcmb4, gamma_isrf, gasgr, nbirthcloud);
            //if (gp->verbose) printf("sol: td=%g tg=%g sol=%g tdold=%g %g %g\n",td, gp->tgas, sol, tdold, tdlo, tdhi);
	    if (sol > 0) tdlo = td;  // heating, so tdust should increase
	    else tdhi = td;
	    if (fabs(tdhi-tdlo) < tol * (tdhi+tdlo) || tdlo > tdhi) break;
	    if (td > gp->tgas || td < tcmb) break;
            //if (gp->verbose) printf("sol: sol=%g nH=%g td=%g tg=%g tdold=%g %g %g isrf=%g\n",sol, gp->nH, td, gp->tgas, tdold, tdlo, tdhi, gamma_isrf);
            if (++iter > maxiter) {
                printf("Crackle: Non-convergence in calculate_dust_temp(), returning tdust=%g (tdold=%g, tdlo=%g tdhi=%g tgas=%g tbc=%g tcmb=%g)\n",td, tdold, tdlo, tdhi, gp->tgas, tbirthcloud, tcmb);
                break;
            }
        }

	gp->tdust = td;
	//assert((gp->tdust >= tdlo) && (gp->tdust <= tdhi));
}

static inline double calculate_dust_temp(double tgas, double nh, double gasgr, double gamma_isrf, double tcmb, double td)
{
        /* Solve for Tdust from ISRF heating, CMB heating, and gas-grain cooling */
        const double tol = 1.e-3;
        double eps=1.e-2, tdold = -1.e20, dsol, slope; // for Newton-Raphson
        //double tdlo=tcmb, tdhi=tgas; // bisection limits
        const double tcmb4 = tcmb * tcmb * tcmb * tcmb;
        int iter = 0, maxiter = 100;
        double sol;

        /* Solve for Tdust via Newton-Raphson */
        if (td<=0.f) td = dust_thermal_balance(td, tgas, tcmb, tcmb4, gamma_isrf, gasgr, nh); // initial guess
        while (fabs(td-tdold) > tol * td) {
            sol = dust_thermal_balance(td, tgas, tcmb, tcmb4, gamma_isrf, gasgr, nh);
            dsol = dust_thermal_balance((1.f+eps)*td, tgas, tcmb, tcmb4, gamma_isrf, gasgr, nh) - sol;
            tdold = td;
	    if (dsol > 0.) td = td - td * eps * sol / dsol;  // Newton-Raphson guess
	    else break;  // converged
            if (td < tcmb) td = tcmb;  // Limit tdust to [tCMB, tgas]
            if (td > tgas) td = tgas;
            if (sol * (sol+dsol) < 0.f) eps *= 0.5f;  // we have passed minimum; reduce eps
	    if (td > tgas || td < tcmb) break;
            if (++iter > maxiter) {
                printf("Crackle: Non-convergence in calculate_dust_temp(), returning tdust=%g (tdold=%g, dsol=%g eps=%g slope=%g)\n",td, tdold, dsol, eps, slope);
                break;
            }
        }

        return td;
}

