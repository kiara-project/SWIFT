
/* Main driver for Crackle chemistry solver, which is a fully C version of Grackle .
 * The code explicitly evolves the internal energy, electron density, and HI density.
 * Other quantities are evolved via a predictor-corrector scheme, if explicit evolution is not accurate enough.
 *
 * Romeel Dave, Feb 2024
 * */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <float.h>
#include "phys_constants.h"
/* Grackle includes */
//#include "grackle.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
/* Crackle-specific includes */
#include "crackle.h"
#include "crackle_aux.h"
#include "crackle_interpolation.h"
#include "crackle_dust.h"

#define MINFRAC 1.e-3
#define CONVERGENCE 0.1  // convergence tolerance, in units of chemistry->accuracy
#define HEATLIM 1000.  // max heating factor in T allowed in a single system timestep
#define COOLLIM 1.e-3  // max cooling factor in T allowed in a single system timestep
#define ITER_HEATLIM 10.  // max heating factor in u allowed in a single iteration
#define JAC_EPS 0.05  // fractional change in u for the finite-difference d(edot)/du; coarse to tolerate tdust solver noise
#define H2_FRAC_MIN_NONISM 1.e-3  // outside the ISM, H2/H mass fraction below which H2 is returned to HI rather than evolved


int crackle_solve_chemistry(grackle_field_data *p, chemistry_data *chemistry, chemistry_data_storage grackle_rates, photo_rate_storage my_uvb_rates, code_units *units, double dt) 
{
	int iter=0;
	double dtit=0., dtcool=0., dtsuppress=1., edot_prev=0.;
	grackle_part_data gp, gp_old, gp_species;
	chemistry_rate_storage my_rates;  // interpolated rates for this field
	interp_struct interpolation;
	crackle_units cunits;

	/* Copy to grackle_part_data */
	copy_grackle_fields_to_part(p, &gp, chemistry);

	check_hydrogen(&gp, 0);

	/* Set up various unit conversions etc */
	set_crackle_units(units, grackle_rates, chemistry->Gamma, &cunits);

	/* Compute some basic properties, save T at start */
	set_rhot(&gp, units, chemistry);
	gp.tstart = gp.tgas;

	/* ISM flag: see whether we should be evolving dust and H2 */
	int ism_flag = (gp.isrf_habing >= 0.);
	if (gp.isrf_habing < 0.) gp.isrf_habing = 0.;

	/* initialize inteprolation of chemistry rate tables for this particle */
	init_temperature_interpolation(&gp, chemistry, &interpolation, &cunits, &grackle_rates, &my_uvb_rates);

	while (dtcool < dt) {
	    //if (gp.tgas < 100.) gp.verbose = 1;
	    //if (gp.id == 3534105) gp.verbose = 1;
	    //if ( (gp.tgas > 1.e6) ) gp.verbose = 1;
	    /* Retain previous iteration particle info */
	    memmove(&gp_old, &gp, sizeof(gp)); 
	    /* Set up cooling/heating rates interpolation */
	    setup_temperature_interpolation(gp.logtem, chemistry, &interpolation);
	    /* Get interpolated chemistry rates for this particle */
	    lookup_chemistry_coeffs(chemistry->primordial_chemistry, &grackle_rates, &my_rates, &interpolation);
	    /* In the first iteration we need edot for the timestep; afterwards we use the effective edot of the previous iteration */
	    if (iter == 0) {
	        compute_edot(&gp, chemistry, &grackle_rates, &my_rates, &my_uvb_rates, &interpolation, units, &cunits, ism_flag);
	    }
	    /* Species rates for the timestep must be recomputed every iteration, else a large initial rate keeps dtit tiny all step */
	    gp.HIdot = compute_HIdot(chemistry->primordial_chemistry, &gp, chemistry, &my_rates, units);
	    gp.dedot = compute_dedot(chemistry->primordial_chemistry, &gp, chemistry, &my_rates, units);
	    /* Set the timestep for this iteration (lowering rates if they are too large) */
	    float accuracy = chemistry->accuracy;
	    //if (ism_flag) accuracy = 0.5f;
	    dtit = compute_iteration_dt(&gp, edot_prev, chemistry, dt, dtcool, accuracy, &dtsuppress);

	    /* If evol is driven by external heating/cooling there is no point iterating because edot won't change much */
	    if (!ism_flag && fabs(gp.edot_ext) > EDOT_EXT_FACTOR * fabs(gp.edot-gp.edot_ext) && dtit < 0.25 * (dt-dtcool)) {
		dtit = dt - dtcool;
	    }

	    /* Evolve H, He species densities */
	    evolve_hydrogen(&gp, &gp, chemistry, &my_rates, dtit);  
	    evolve_helium(&gp, &gp, chemistry, &my_rates, dtit); 
	    if (chemistry->primordial_chemistry >= 2) evolve_H2(&gp, ism_flag, chemistry, &my_rates, &cunits, dtit);  
	    evolve_elements(&gp, &gp_old, chemistry);

	    /* Save step-averaged species densities so we can redo iteration at higher accuracy */
	    evolve_pred_corr(&gp, &gp_old, chemistry);
	    save_species_data(&gp, &gp_species);
	    memmove(&gp, &gp_old, sizeof(gp)); 
	    save_species_data(&gp_species, &gp);

	    /* Energy update with step-averaged species: 2nd-order L-stable Rosenbrock (ROS2) step, stable even when dtit is much
	     * longer than the thermal relaxation time (e.g. near thermal equilibrium, where |edot| << individual heating/cooling terms) */
	    evolve_energy_ros2(&gp, gp_old.internal_energy, dtit, chemistry, &grackle_rates, &my_rates, &my_uvb_rates, &interpolation, units, &cunits, ism_flag);

	    /* Species rates at the mid-step temperature for the second species half-step */
	    grackle_part_data gp_mid = gp;
	    gp_mid.internal_energy = 0.5 * (gp_old.internal_energy + gp.internal_energy);
	    set_rhot(&gp_mid, units, chemistry);
	    setup_temperature_interpolation(gp_mid.logtem, chemistry, &interpolation);
	    lookup_chemistry_coeffs(chemistry->primordial_chemistry, &grackle_rates, &my_rates, &interpolation);

	    /* Do iteration again using step-averaged species densities; only need 1/2 iteration for species */
	    evolve_hydrogen(&gp, &gp, chemistry, &my_rates, 0.5 * dtit);  
	    evolve_helium(&gp, &gp, chemistry, &my_rates, 0.5 * dtit); 
	    if (chemistry->primordial_chemistry >= 2) evolve_H2(&gp, ism_flag, chemistry, &my_rates, &cunits, 0.5 * dtit);  
	    if (chemistry->use_dust_evol) evolve_dust(&gp, chemistry, units, ism_flag, dtit); 
	    evolve_elements(&gp, &gp_old, chemistry);

	    /* Done with iteration */
	    set_rhot(&gp, units, chemistry);
	    dtcool += dtit;
	    edot_prev = gp.edot;
	    iter ++;

	    if (gp.verbose) printf("iter: i=%d dt=%g fdt=%g id=%d nh=%g e=%g de=%g HI=%g HII=%g H2I=%g HeI=%g HeII=%g dust=%g T=%g Td=%g\n", iter, dtit, dtcool/dt, gp.id, gp.rhoH * units->density_units / mh, gp.internal_energy, gp.e_density, gp.HI_density/gp.density, gp.HII_density/gp.density, gp.H2I_density/gp.density, gp.HeI_density/gp.density, gp.HeII_density/gp.density, gp.dust_density/gp.density, gp.tgas, gp.tdust);
	    if (gp.verbose) printf("rates: i=%d fdt=%g id=%d edot=%g dedot=%g HIdot=%g T=%g n=%g\n",iter, dtit/dt, gp.id, gp.edot, gp.dedot, gp.HIdot, gp.tgas, gp.density*units->density_units/mh );
	   
	    gp.verbose = 0; 

	    /* If we've reached temp floor/ceiling and are still cooling/heating, then we're done */
	    if (apply_temperature_bounds(&gp, chemistry, fmax(gp.temperature_floor, COOLLIM * gp.tstart), HEATLIM * gp.tstart)) break;

	    /* If key quantities are converged then we're done. The change in this iteration is projected over the
	     * remaining time, since a small change over a small dtit does not mean the system is in equilibrium. */
	    const double fproj = (dt - dtcool) / dtit;
	    const double conv_tol = CONVERGENCE * chemistry->accuracy;
	    if (fabs(gp_old.internal_energy - gp.internal_energy) * fproj < conv_tol * gp_old.internal_energy
		&& fabs(gp_old.HI_density - gp.HI_density) * fproj < conv_tol * fmax(gp_old.HI_density, MINFRAC*gp_old.density)
		&& fabs(gp_old.e_density - gp.e_density) * fproj < conv_tol * fmax(gp_old.e_density, MINFRAC*gp_old.density) ) break;

	    /* Too many iterations */
	    if (iter > chemistry->max_iterations) break;
	}

	/* If there is time left over, need to form/destroy H2 and dust over the full timestep */
	if (dtcool < dt && chemistry->primordial_chemistry >= 2) {
	    iter = chemistry->max_iterations;
	    evolve_hydrogen(&gp, &gp, chemistry, &my_rates, dt - dtcool);  
	    evolve_H2(&gp, ism_flag, chemistry, &my_rates, &cunits, dt - dtcool); 
	    evolve_elements(&gp, &gp_old, chemistry);
	    if (chemistry->use_dust_evol) evolve_dust(&gp, chemistry, units, ism_flag, dt - dtcool);
	}

	/* Copy from grackle_part_data */
	copy_grackle_fields_from_part(p, &gp, chemistry);
	check_hydrogen(&gp, 99);

	return 1;
}

double compute_iteration_dt(grackle_part_data *gp, double edot_prev, chemistry_data *chemistry, double dt, double dtcool, float accuracy, double *dtsuppress)
{
	/* Calculate timesteps for each tracked quantity */
	/* edot_end catches accelerating cooling/heating, which the effective edot of the last iteration lags */
	const double dt_e = fabs( gp->internal_energy * gp->density / (fmax(fabs(gp->edot), fabs(gp->edot_end)) + tiny) );
	const double dt_de = fabs( fmax(gp->e_density, MINFRAC*gp->density) / (gp->dedot + tiny) );
	const double dt_HI = fabs( fmax(gp->HI_density, MINFRAC*gp->density) / (gp->HIdot + tiny) );
	double dt_H2 = FLT_MAX;
	if (chemistry->primordial_chemistry >= 2) {
	    dt_H2 = fabs( fmax(gp->H2I_density, MINFRAC*gp->density) / (gp->H2Idot + tiny) );
	}

	/* If thermal equilbirium is passed (edot changed sign since last iteration), attenuate the timestep to help it converge */
	if (gp->edot * edot_prev < 0.) *dtsuppress *= 0.5;
	else *dtsuppress *= (1. + chemistry->accuracy);  // asymmetric wrt the reduction to avoid oscillations
	if (*dtsuppress > 1.) *dtsuppress = 1.;  // only ever reduce dtit below the accuracy-limited value

	/* Compute the timestep for this iteration */
	double dtit = chemistry->accuracy * *dtsuppress * fmin(fmin(dt_de, dt_HI), dt_e);
	if (chemistry->primordial_chemistry >= 2) {
	    dtit = fmin(dtit, chemistry->accuracy * dt_H2);
	}

	/* Limit dtit if we will reach end of timestep this iteration */
	if (dtit > dt - dtcool) {
	    dtit = dt - dtcool + tiny;
	}

	return dtit;
}

/* du/dt = edot/density at internal energy u, holding species fixed; updates the derived quantities and rates in g, rates, interp */
static double crackle_dudt(grackle_part_data *g, double u, chemistry_data *chemistry, const chemistry_data_storage *grackle_rates, chemistry_rate_storage *rates, const photo_rate_storage *my_uvb_rates, interp_struct *interp, code_units *units, const crackle_units *cunits, int ism_flag)
{
	g->internal_energy = u;
	set_rhot(g, units, chemistry);
	setup_temperature_interpolation(g->logtem, chemistry, interp);
	lookup_chemistry_coeffs(chemistry->primordial_chemistry, grackle_rates, rates, interp);
	compute_edot(g, chemistry, grackle_rates, rates, my_uvb_rates, interp, units, cunits, ism_flag);
	return g->edot / g->density;
}

void evolve_energy_ros2(grackle_part_data *gp, double u0, double dtit, chemistry_data *chemistry, const chemistry_data_storage *grackle_rates, chemistry_rate_storage *my_rates, const photo_rate_storage *my_uvb_rates, interp_struct *interpolation, code_units *units, const crackle_units *cunits, int ism_flag)
{
	/* ROS2 (Verwer et al. 1999): (1-g h J) k1 = f(u0); (1-g h J) k2 = f(u0 + h k1) - 2 k1; u1 = u0 + h (1.5 k1 + 0.5 k2).
	 * Second order for any approximation to J, so a coarse finite-difference J is fine. */
	const double gamma_ros = 1. + 1./sqrt(2.);

	/* f(u0); evaluated on gp itself so tdust, shielding factors, edot_ext etc. are retained */
	const double f0 = crackle_dudt(gp, u0, chemistry, grackle_rates, my_rates, my_uvb_rates, interpolation, units, cunits, ism_flag);

	/* Jacobian by finite difference, using scratch copies. Only the damping part (J<0) is used; where J>0 (thermally unstable,
	 * e.g. cooling decreasing with T) the step reduces to explicit 2nd-order Heun. */
	/* The sub-evaluations change only u, so the self-shielding (weakly T-dependent) is held at its f(u0) value, like the species */
	grackle_part_data g_tmp = *gp;
	g_tmp.freeze_shielding = 1;
	chemistry_rate_storage rates_tmp = *my_rates;
	interp_struct interp_tmp = *interpolation;
	const double fJ = crackle_dudt(&g_tmp, (1. + JAC_EPS) * u0, chemistry, grackle_rates, &rates_tmp, my_uvb_rates, &interp_tmp, units, cunits, ism_flag);
	const double jac = fmin((fJ - f0) / (JAC_EPS * u0), 0.);
	const double w = 1. / (1. - gamma_ros * dtit * jac);

	/* Stage 2, with the predicted energy kept within the per-iteration limits */
	const double k1 = w * f0;
	double u_pred = u0 + dtit * k1;
	u_pred = fmin(fmax(u_pred, fmax(u0 / ITER_HEATLIM, gp->u_cmb)), ITER_HEATLIM * u0);
	g_tmp = *gp;
	g_tmp.freeze_shielding = 1;
	rates_tmp = *my_rates;
	interp_tmp = *interpolation;
	const double f1 = crackle_dudt(&g_tmp, u_pred, chemistry, grackle_rates, &rates_tmp, my_uvb_rates, &interp_tmp, units, cunits, ism_flag);
	const double k2 = w * (f1 - 2. * k1);
	gp->edot_end = gp->density * f1;

	/* Update u from u0 with the effective edot, applying the usual per-iteration limits */
	gp->internal_energy = u0;
	gp->edot = gp->density * (1.5 * k1 + 0.5 * k2);
	evolve_internal_energy(gp, chemistry, dtit);
	if (gp->verbose) printf("ros2: id=%d u0=%g f0=%g f1=%g jac*dt=%g u1/u0=%g\n", gp->id, u0, f0, f1, jac*dtit, gp->internal_energy/u0);
}

void evolve_internal_energy(grackle_part_data *gp, chemistry_data *chemistry, double dtit) 
{
	const double u_prev = gp->internal_energy;
	gp->internal_energy += gp->edot / gp->density * dtit;
	/* Limits -- don't let u change by more than accuracy level in a single iteration */
	if (gp->internal_energy < (1.-chemistry->accuracy) * u_prev) gp->internal_energy = (1.-chemistry->accuracy) * u_prev;
	if (gp->internal_energy > ITER_HEATLIM * u_prev) gp->internal_energy = ITER_HEATLIM * u_prev;
	if (gp->internal_energy < gp->u_cmb) gp->internal_energy = gp->u_cmb;

	return;
}

void compute_edot(grackle_part_data *gp, chemistry_data *chemistry, const chemistry_data_storage *grackle_rates, chemistry_rate_storage *my_rates, const photo_rate_storage *my_uvb_rates, interp_struct *interpolation, code_units *units, const crackle_units *cunits, int ism_flag)
{
	double edot_prim = 0., edot_h2 = 0., edot_gasgr = 0., edot_uvb = 0., edot_pe = 0., edot_edust = 0., edot_comp = 0., edot_rt = 0., edot_h2heat = 0., edot_ext = 0., edot_metal = 0.;

	/* Interpolate primordial excitation, ionization, and recombination rates for given T */
	primordial_cooling_rates(interpolation, grackle_rates, my_rates);

	gp->edot = 0.;  // initialize for this step
	gp->edot_ext = 0.;  // must be reset each call, else it accumulates over iterations

	/* Primordial cooling */
	edot_prim += -my_rates->ceHI*gp->HI_density*gp->e_density 
		- 0.25*my_rates->ceHeI*gp->HeII_density*gp->e_density*gp->e_density*cunits->dom 
		- 0.25*my_rates->ceHeII*gp->HeII_density*gp->e_density;
	edot_prim += -my_rates->ciHI*gp->HI_density*gp->e_density
		- 0.25*my_rates->ciHeI*gp->HeI_density*gp->e_density
		- 0.25*my_rates->ciHeII*gp->HeII_density*gp->e_density
		- 0.25*my_rates->ciHeIS*gp->HeII_density*cunits->dom*gp->e_density*gp->e_density;
	edot_prim += -my_rates->reHII* gp->HII_density*gp->e_density
		- 0.25*my_rates->reHeII1*gp->HeII_density*gp->e_density
		- 0.25*my_rates->reHeII2*gp->HeII_density*gp->e_density
		- 0.25*my_rates->reHeIII*gp->HeIII_density*gp->e_density;
	edot_prim -= my_rates->brem*(gp->HII_density+0.25*gp->HeII_density+gp->HeIII_density)*gp->e_density;
	gp->edot += edot_prim;

	/* Calculate self-shielding */
	if (!gp->freeze_shielding) gp->f_shield = compute_self_shielded_rates(gp, chemistry, my_rates, my_uvb_rates, cunits, interpolation);
	if (ism_flag == 0) gp->f_shield = 1.f;

	/* Add H2 cooling (Glover & Abel 2008); fudge at extreme density NOT included */
	if (chemistry->primordial_chemistry >= 2) {
	    H2_cooling_rates(interpolation, grackle_rates, my_rates);
	    const double galdl = my_rates->GAHI * gp->HI_density + 0.5 * my_rates->GAH2 * gp->H2I_density + 0.25 * my_rates->GAHe * gp->HeI_density + my_rates->GAHp * gp->HII_density + my_rates->GAel * gp->e_density;
	    const double gphdl1 = my_rates->H2LTE * cunits->dom_inv;
	    edot_h2 -= gp->H2I_density * my_rates->H2LTE / ((1.+gphdl1/galdl) * 2.*cunits->dom);
	    gp->edot += edot_h2;

	    /* Dust-related cooling/heating processes */
	    if (crackle_h2_on_dust(chemistry)) {
	        /* compute dust temperature */
	        double gasgr = interpolate_rates(grackle_rates->gas_grain, interpolation->binfrac, interpolation->index);
	        double gasgr_tdust = chemistry->local_dust_to_gas_ratio * gasgr * cunits->coolunit / mh;
	        /* Get dust rates */
	        //gp->tdust = calculate_dust_temp(gp->tgas, gp->nH, gasgr_tdust, grackle_rates->gamma_isrf * gp->isrf_habing, cunits->compton2, gp->tdust); // calc_tdust_1d_g()
	        /* ISRF dust heating is not attenuated by f_shield: that is the H2 Lyman-Werner self-shielding factor, which does not apply to continuum absorption by grains */
	        calculate_tdust_bisect(gp, gasgr_tdust, grackle_rates->gamma_isrf * gp->isrf_habing, cunits->compton2);
		//if(gp->tdust < 2*cunits->compton2) printf("DUST td=%g gasgr=%g tcmb=%g isrf=%g\n",gp->tdust, gasgr_tdust, cunits->compton2, gp->isrf_habing);
	        /* Gas-dust grain heat transfer rate */
	        dust_species_rates(gp->tdust, gp->dust2gas, chemistry, grackle_rates, my_rates, interpolation);
	        if (gp->dust2gas > tiny) {
	            edot_gasgr -= gasgr * (gp->tgas - gp->tdust) * gp->dust2gas * gp->rhoH * gp->rhoH;
	            gp->edot += edot_gasgr;
	        }
	    }
	}

	/* UVB heating */
	if (chemistry->UVbackground && gp->grid_end >= 0) {
	    edot_uvb += cunits->dom_inv * (my_uvb_rates->piHI * gp->fSShHI * gp->HI_density + 0.25 * (my_uvb_rates->piHeI * gp->fSShHeI * gp->HeI_density + my_uvb_rates->piHeII * gp->fSShHeII * gp->HeII_density));
	    gp->edot += edot_uvb;
	}

	/* Photo-electric heating by UV-irradiated dust */
	if (chemistry->use_isrf_field && gp->isrf_habing > tiny && gp->dust2gas > tiny) {
	    double gammah_eff = my_rates->gammah * 0.05 * gp->isrf_habing;  // Wolfire+1995 eq 1
	    if (gp->tgas > 2.e4) gammah_eff = 0.;
	    edot_pe += gammah_eff * gp->rhoH * cunits->dom_inv * gp->dust2gas / chemistry->local_dust_to_gas_ratio;
	    gp->edot += edot_pe;
	    //printf("pe %g %g %g\n",gp->edot,gp->dust2gas,gammah_eff);
	}

	/* Electron recombination onto dust grains */
	if ((chemistry->use_dust_evol || chemistry->dust_chemistry) && gp->metallicity > tiny && gp->e_density > tiny && gp->isrf_habing > tiny) {
	    const double regr = interpolate_rates(grackle_rates->regr, interpolation->binfrac, interpolation->index);
	    const double grbeta = 0.74 * pow(gp->tgas, -0.068);
	    edot_edust -= regr * pow(gp->isrf_habing * cunits->dom_inv / gp->e_density, grbeta) * gp->e_density * gp->rhoH * gp->dust2gas / chemistry->local_dust_to_gas_ratio;
	    gp->edot += edot_edust;
	}

	/* CMB Compton cooling and X-ray compton heating */
	edot_comp -= cunits->compton1 * (gp->tgas-cunits->compton2) * gp->e_density * cunits->dom_inv;
	edot_comp -= my_uvb_rates->comp_xray * (gp->tgas-my_uvb_rates->temp_xray) * gp->e_density * cunits->dom_inv;
	gp->edot += edot_comp;

	/* Photoheating from radiative transfer */
	if (chemistry->use_radiative_transfer && gp->RT_heating_rate > FLT_MIN) {
	    edot_rt += gp->RT_heating_rate * gp->fSShHI * gp->HI_density * cunits->dom_inv * cunits->coolunit_inv;
	    gp->edot += edot_rt;
	    gp->edot_ext += edot_rt;
	}

	/* H2 heating */
	if (chemistry->primordial_chemistry >= 2) {
	    double h2heatfac = 1. + my_rates->n_cr_n / (cunits->dom * (my_rates->n_cr_d1 * gp->HI_density + 0.5 * my_rates->n_cr_d2 * gp->H2I_density));
	    h2heatfac = 1./h2heatfac;
	    double H2delta = 4.48 * gp->HI_density * (my_rates->k22 * gp->HI_density * gp->HI_density - 0.5 * my_rates->k13 * gp->H2I_density);
	    if (H2delta > tiny) H2delta *= h2heatfac;
	    if (crackle_h2_on_dust(chemistry) && gp->dust2gas > tiny) {
		H2delta += my_rates->h2dust * gp->HI_density * gp->rhoH * (0.2 + 4.2 * h2heatfac);
	    }
	    edot_h2heat += cunits->chunit * H2delta;
	    gp->edot += edot_h2heat;
	}

	/* Cloudy-based metal cooling/heating */
	if (chemistry->metal_cooling && gp->metallicity > tiny) {
	    edot_metal += gp->rhoH * gp->rhoH * gp->metallicity * cloudy_metal_cooling(gp->logtem, chemistry, interpolation, grackle_rates, cunits);  // cool1d_cloudy_g()
	    gp->edot += edot_metal;
	}

	/* User-specified external heating terms */
	if (chemistry->use_volumetric_heating_rate && gp->volumetric_heating_rate > tiny) {
	    edot_ext += gp->volumetric_heating_rate * cunits->dom_inv * cunits->dom_inv * cunits->coolunit_inv;
	}
	if (chemistry->use_specific_heating_rate && gp->specific_heating_rate > tiny) {
	    edot_ext += gp->specific_heating_rate * gp->density * mh * cunits->dom_inv * cunits->coolunit_inv;
	}
	gp->edot += edot_ext;
	gp->edot_ext += edot_ext;

	//gp->edot = edot_prim + edot/_h2 + edot_gasgr + edot_uvb + edot_pe + edot_edust + edot_comp + edot_rt + edot_h2heat + edot_ext + edot_metal;
	if (gp->verbose) {
	    printf("edot: id=%d nH=%g isrf=%g tgas=%g told=%g edot=%g pr=%g h2=%g gr=%g uvb=%g pe=%g ed=%g co=%g rt=%g %g h2h=%g ext=%g met=%g\n", gp->id, gp->nH, gp->isrf_habing, gp->tgas, gp->tstart, gp->edot, edot_prim , edot_h2 , edot_gasgr , edot_uvb, edot_pe , edot_edust , edot_comp , edot_rt , gp->RT_heating_rate, edot_h2heat , edot_ext , edot_metal); 
	    fflush(stdout);
	}

	assert(gp->edot==gp->edot && "CRACKLE ERROR: gp->edot is nan.");  // check for NaN

	return;
}

double compute_dedot(int chemistry_flag, const grackle_part_data *gp, chemistry_data *chemistry, const chemistry_rate_storage *my_rates, code_units *units) 
{ 
	double dedot=0.;

	/* Recombination */
	dedot = my_rates->k1 * gp->HI_density * gp->e_density 
		+ 0.25 * my_rates->k3 * gp->HeI_density * gp->e_density 
		+ 0.25 * my_rates->k5 * gp->HeII_density * gp->e_density
		- my_rates->k2 * gp->HII_density * gp->e_density
		- 0.25 * my_rates->k4 * gp->HeII_density * gp->e_density
		- 0.25 * my_rates->k6 * gp->HeIII_density * gp->e_density
		+ my_rates->k57 * gp->HI_density * gp->HI_density 
		+ 0.25 * my_rates->k58 * gp->HI_density * gp->HeI_density;

	if (chemistry->UVbackground > 0) {
	    dedot += my_rates->k24shield * gp->HI_density
		+ 0.25 * my_rates->k25shield * gp->HeII_density
		+ 0.25 * my_rates->k26shield * gp->HeI_density;
	}
	/* RT photoionization */
	if (chemistry->use_radiative_transfer && gp->RT_heating_rate > FLT_MIN) {
	    dedot += gp->RT_HI_ionization_rate * gp->fSShHI * gp->HI_density + 0.25 * (gp->RT_HeI_ionization_rate * gp->fSShHeI * gp->HeI_density + gp->RT_HeII_ionization_rate * gp->fSShHeII * gp->HeII_density);
	}

	/* Molecular */
	if (chemistry_flag > 1) {
	    dedot += my_rates->k8 * gp->HI_density * gp->HM_density
		+ my_rates->k15 * gp->HI_density * gp->HM_density
		+ my_rates->k17 * gp->HII_density * gp->HM_density
		+ my_rates->k14 * gp->e_density * gp->HM_density
	    	- my_rates->k7 * gp->HI_density * gp->e_density 
		- 0.5 * my_rates->k18 * gp->H2II_density * gp->e_density;
	}

	return dedot;
}

double compute_HIdot(int chemistry_flag, const grackle_part_data *gp, chemistry_data *chemistry, const chemistry_rate_storage *my_rates, code_units *units) 
{ 
	double HIdot=0.;

	/* Recombination */
	HIdot = -my_rates->k1* gp->HI_density * gp->e_density 
		+ my_rates->k2* gp->HII_density * gp->e_density 
	/* Collisional ionization */
		- my_rates->k57 * gp->HI_density * gp->HI_density 
		- 0.25 * my_rates->k58 * gp->HI_density * gp->HeI_density 
	/* UVB photo-ionization (no self-shielding -- assumes self-shielding built into ionization table if desired */
		- my_rates->k24shield * gp->HI_density;
	/* RT photoionization */
	if (chemistry->use_radiative_transfer) {
	    HIdot -= gp->RT_HI_ionization_rate * gp->fSShHI * gp->HI_density;
	}

	/* Molecular */
	if (chemistry_flag >= 2) {
	    HIdot += -my_rates->k7 * gp->e_density * gp->HI_density 
		    - my_rates->k8 * gp->HM_density * gp->HI_density 
		    - my_rates->k9 * gp->HII_density * gp->HI_density 
		    - 0.5 *my_rates->k10 * gp->H2II_density * gp->HI_density 
		    - 2.f *my_rates->k22 * gp->HI_density * gp->HI_density * gp->HI_density
		    + 0.5 *my_rates->k11 * gp->HII_density * gp->H2I_density 
		    + my_rates->k12 * gp->e_density * gp->H2I_density 
		    + my_rates->k13 * gp->H2I_density * gp->HI_density 
		    + my_rates->k14 * gp->e_density * gp->HM_density 
		    + my_rates->k15 * gp->HM_density * gp->HI_density
		    + 2.f *my_rates->k16 * gp->HII_density * gp->HM_density 
		    + my_rates->k18 * gp->e_density * gp->H2II_density
		    + 0.5 * my_rates->k19 * gp->HM_density * gp->H2II_density;
	    HIdot += my_rates->k31shield * gp->H2I_density;  // UVB and/or RT H2 dissociation
	    if (crackle_h2_on_dust(chemistry) && gp->dust2gas > tiny) {
		HIdot -= 2.f *my_rates->h2dust * gp->HI_density * gp->rhoH;
	    }
	}

	return HIdot;
}


void evolve_helium(grackle_part_data *p, grackle_part_data *gp_old, chemistry_data *chemistry, const chemistry_rate_storage *my_rates, double dtit) /* from step_rate_g() */
{ 
	double scoef, acoef;

	/* HeI */
	scoef = my_rates->k4 * p->HeII_density * p->e_density;
	acoef = my_rates->k3 * p->e_density;
	if (chemistry->UVbackground > 0) {
	    acoef += my_rates->k26shield;
	}
	if (chemistry->use_radiative_transfer) acoef += p->RT_HeI_ionization_rate * p->fSShHeI;
	double HeIp = (scoef * dtit + p->HeI_density) / (1.f + acoef * dtit);
	if (HeIp < 0. ) HeIp = 0.;
	p->delta_HeI = HeIp - p->HeI_density;

	/* HeII */
	scoef = my_rates->k3 * HeIp * p->e_density +
		my_rates->k6 * p->HeIII_density * p->e_density;
	acoef = (my_rates->k4 + my_rates->k5) * p->e_density;
	if (chemistry->UVbackground > 0) {
	    scoef += my_rates->k26shield * HeIp;
	    acoef += my_rates->k25shield;
	}
	if (chemistry->use_radiative_transfer) {
	    scoef += p->RT_HeI_ionization_rate * p->fSShHeI * HeIp;
	    acoef += p->RT_HeII_ionization_rate * p->fSShHeII;
	}
	double HeIIp = (scoef * dtit + p->HeII_density) / (1.f + acoef * dtit);
	if (HeIIp < 0. ) HeIIp = 0.;
	p->delta_HeII = HeIIp - p->HeII_density;

	/* HeIII */
	scoef = my_rates->k5 * HeIIp * p->e_density;
	acoef = my_rates->k6 * p->e_density;
	if (chemistry->UVbackground > 0) {
	    scoef += my_rates->k25shield * HeIIp;
	}
	if (chemistry->use_radiative_transfer) scoef += p->RT_HeII_ionization_rate * p->fSShHeII * HeIIp;
	p->delta_HeIII = (scoef * dtit + p->HeIII_density) / (1.f + acoef * dtit) - p->HeIII_density; 
	if (p->delta_HeIII + p->HeIII_density < 0.f ) p->delta_HeIII = -p->HeIII_density;
	/*if (p->delta_HeII != p->delta_HeII || p->delta_HeIII != p->delta_HeIII) {
	    fprintf(stdout, "TROUBLE: dHeII=%g dHeIII=%g scoef=%g acoef=%g HeIIp=%g e=%g k3=%g k4=%g k5=%g k6=%g k25=%g k25shielf=%g k26=%g\n", p->delta_HeII, p->delta_HeIII, scoef, acoef, HeIIp, p->e_density, my_rates->k3, my_rates->k4, my_rates->k5, my_rates->k6, my_rates->k25, my_rates->k25shield, my_rates->k26);
	    fflush(stdout);
	    assert(p->delta_HeII == p->delta_HeII);
	    assert(p->delta_HeIII == p->delta_HeIII);
	}*/

	return;
}

void evolve_hydrogen(grackle_part_data *p, grackle_part_data *gp_old, chemistry_data *chemistry, const chemistry_rate_storage *my_rates, double dtit) /* from step_rate_g() */
{ 
	double scoef, acoef;

	/* HI */
	scoef = my_rates->k2 * p->HII_density * p->e_density;
	if (chemistry->primordial_chemistry >= 2) {
	    scoef += my_rates->k13 * p->HI_density * p->H2I_density
     	          + 0.5 * my_rates->k11 * p->HII_density * p->H2I_density
     	          + my_rates->k12 * p->e_density * p->H2I_density
     	          + my_rates->k14 * p->HM_density * p->e_density
     	          + my_rates->k15 * p->HM_density * p->HI_density
     	          + 2. * my_rates->k16 * p->HM_density * p->HII_density
     	          + my_rates->k18 * p->H2II_density * p->e_density
     	          + 0.5 * my_rates->k19 * p->H2II_density * p->HM_density;
	    scoef += my_rates->k31shield * p->H2I_density;  // UVB and/or RT H2 dissociation
	}
	acoef = my_rates->k1 * p->e_density
	      +	my_rates->k57 * p->HI_density
	      +	0.25 * my_rates->k58 * p->HeI_density;

	if (chemistry->UVbackground > 0) {
     	      acoef += my_rates->k24shield;
	}
	if (chemistry->use_radiative_transfer) {
	      acoef += p->RT_HI_ionization_rate * p->fSShHI;
	}
	if (chemistry->primordial_chemistry >= 2) {
	    acoef += my_rates->k7 * p->e_density
	          + my_rates->k8 * p->HM_density
	          + my_rates->k9 * p->HII_density
	          + 0.5 * my_rates->k10 * p->H2II_density
	          + 2. * my_rates->k22 * p->HI_density * p->HI_density;
	}
	if (crackle_h2_on_dust(chemistry) && p->dust2gas > tiny) {
	    acoef += 2. * my_rates->h2dust * p->rhoH;
	}
	double HIp = (scoef * dtit + p->HI_density) / (1.f + acoef * dtit);
	if (HIp < 0. ) HIp = 0.;
	p->delta_HI = HIp - p->HI_density;

	/* HII */
	scoef = my_rates->k1 * p->HI_density * p->e_density
	      +	0.5 * my_rates->k10 * p->H2II_density * p->HI_density
	      +	my_rates->k57 * p->HI_density * p->HI_density
	      +	0.25 * my_rates->k58 * p->HeI_density * p->HI_density;
	acoef = my_rates->k2 * p->e_density
	      +	my_rates->k9 * p->HI_density
	      +	0.5 * my_rates->k11 * p->H2I_density
	      +	my_rates->k16 * p->HM_density
	      +	my_rates->k17 * p->HM_density;

	if (chemistry->UVbackground > 0) {
	    scoef += my_rates->k24shield * p->HI_density;
	}
	if (chemistry->use_radiative_transfer) {
	    scoef += p->RT_HI_ionization_rate * p->fSShHI * p->HI_density;
	}
	double HIIp = (scoef * dtit + p->HII_density) / (1.f + acoef * dtit);
	if (HIIp < 0. ) HIIp = 0.;
	p->delta_HII = HIIp - p->HII_density;

	/* electrons */
	scoef = my_rates->k8 * p->HM_density * p->HI_density
	      + my_rates->k15 * p->HM_density * p->HI_density
	      + my_rates->k17 * p->HM_density * p->HII_density
	      +	my_rates->k10 * p->H2II_density * p->HI_density
	      +	my_rates->k57 * p->HI_density * p->HI_density
	      +	0.25 * my_rates->k58 * p->HeI_density * p->HI_density;

	if (chemistry->UVbackground > 0) {
	    scoef += my_rates->k24shield * HIp
     	          + 0.25 * my_rates->k25shield * (p->HeII_density+p->delta_HeII)
     	          + 0.25 * my_rates->k26shield * (p->HeI_density+p->delta_HeI);
	}
	if (chemistry->use_radiative_transfer) {
	    scoef += p->RT_HI_ionization_rate * p->fSShHI*  HIp
	          + 0.25 * p->RT_HeI_ionization_rate * p->fSShHeI * (p->HeI_density+p->delta_HeI)
	          + 0.25 * p->RT_HeII_ionization_rate * p->fSShHeII * (p->HeII_density+p->delta_HeII);
	}
	acoef = my_rates->k1 * p->HI_density
	      - my_rates->k2 * p->HII_density
	      +	0.25 * my_rates->k3 * p->HeI_density
	      -	0.25 * my_rates->k6 * p->HeIII_density
	      +	0.25 * my_rates->k5 * p->HeII_density
	      -	0.25 * my_rates->k4 * p->HeII_density
	      +	my_rates->k14 * p->HM_density
	      -	my_rates->k7 * p->HI_density
	      -	0.5 * my_rates->k18 * p->H2II_density;
	acoef *= -1.;
	p->delta_e = (scoef * dtit + p->e_density) / (1.f + acoef * dtit) - p->e_density;

	return;
}

void evolve_H2(grackle_part_data *p, int ism_flag, chemistry_data *chemistry, const chemistry_rate_storage *my_rates, const crackle_units *cunits, double dtit) /* from step_rate_g() */
{ 
	double scoef, acoef;

	/* Outside the ISM, H2 is evolved only while above a small fraction (e.g. H2 carried out of the ISM);
	 * below this it is returned to atomic H, so H2 does not form from scratch outside the ISM (for speed) */
	const double rhoH2 = p->H2I_density + p->H2II_density;
	if (ism_flag == 0 && rhoH2 < H2_FRAC_MIN_NONISM * p->rhoH) {
	    if (rhoH2 > 0.f || p->HM_density > 0.f) {
	        p->HI_density += p->H2I_density + p->HM_density;
	        p->HII_density += p->H2II_density;
	        p->H2I_density = p->H2II_density = p->HM_density = 0.f;
	        p->H2Idot = 0.f;
	        compute_electron_density(p);
	    }
	    /* clear any changes from a previous iteration when H2 was still being evolved */
	    p->delta_H2I = p->delta_H2II = p->delta_HM = 0.f;
	    return;
	}

	/* H2 */
	scoef = 2.f * my_rates->k8 * p->HM_density * p->HI_density +
		my_rates->k10 * p->H2II_density * p->HI_density +
		my_rates->k19 * p->H2II_density * p->HM_density +
		2.f * my_rates->k22 * p->HI_density * p->HI_density * p->HI_density;
	acoef = my_rates->k13 * p->HI_density + 
		my_rates->k11 * p->HII_density +
		my_rates->k12 * p->e_density;
	acoef += my_rates->k31shield;  // UVB and/or RT H2 dissociation
	if (chemistry->UVbackground > 0) {
	    acoef += my_rates->k29shield;
	}
	if (crackle_h2_on_dust(chemistry) && p->dust2gas > tiny) {
	    scoef += 2.f * my_rates->h2dust * p->HI_density * p->rhoH;
	}
	p->delta_H2I = (scoef * dtit + p->H2I_density) / (1.f + acoef * dtit) - p->H2I_density;
	if (p->delta_H2I < -p->H2I_density) p->delta_H2I = -p->H2I_density;

	/* H- */
	scoef = my_rates->k7 * p->e_density * p->HI_density;
	acoef = (my_rates->k8 + my_rates->k15) * p->HI_density + 
		(my_rates->k16 + my_rates->k17) * p->HII_density + 
		my_rates->k14 * p->e_density +
		0.5 * my_rates->k19 * p->H2II_density;
	if (chemistry->UVbackground > 0) {
	    acoef += my_rates->k27;
	}
	p->delta_HM = (scoef * dtit + p->HM_density) / (1.f + acoef * dtit) - p->HM_density;

	/* H2+ */
	/* For some reason H2+ uses the updated densities, so let's calculate those */
	double HIp, HIIp, H2Ip, HMp, dep;
	HIp = p->HI_density + p->delta_HI;
	HIIp = p->HII_density + p->delta_HII;
	H2Ip = p->H2I_density + p->delta_H2I;
	HMp = p->HM_density + p->delta_HM;
	dep = p->e_density + p->delta_e;
	/* Equilibrium H2+ = formation / destruction rate. The destruction rate can be exactly zero, e.g. fully molecular
	 * neutral gas (HI = e = H- = 0) above the UVB redshift range; then H2+ has no reactions and is left unchanged. */
	double H2II_form = my_rates->k9 * HIp * HIIp + 0.5f * my_rates->k11 * H2Ip * HIIp + my_rates->k17 * HMp * HIIp;
	double H2II_destr = my_rates->k10 * HIp + my_rates->k18 * dep + my_rates->k19 * HMp;
	if (chemistry->UVbackground > 0) {
	    H2II_form += my_rates->k29shield * H2Ip;
	    H2II_destr += my_rates->k28shield + my_rates->k30shield;
	}
	if (H2II_destr > 0.) p->delta_H2II = 2.f * H2II_form / H2II_destr - p->H2II_density;
	else p->delta_H2II = 0.;

	/*if (p->delta_H2I != p->delta_H2I || p->delta_H2II != p->delta_H2II) {
	    fprintf(stdout, "TROUBLE: dH2I=%g dH2II=%g scoef=%g acoef=%g H2Ip=%g e=%g k3=%g k4=%g k5=%g k6=%g k25=%g k25shielf=%g k26=%g\n", p->delta_H2I, p->delta_H2II, scoef, acoef, H2Ip, p->e_density, my_rates->k3, my_rates->k4, my_rates->k5, my_rates->k6, my_rates->k25, my_rates->k25shield, my_rates->k26);
	    fflush(stdout);
	    assert(p->delta_H2I == p->delta_H2I);
	    assert(p->delta_H2II == p->delta_H2II);
	}*/

	p->H2Idot = p->delta_H2I / dtit;

	return;
}

void evolve_elements(grackle_part_data *gp, grackle_part_data *gp_old, chemistry_data *chemistry)
{
	check_hydrogen(gp, 1);
	/* Limit change for H species */
	if (gp->HI_density + gp->delta_HI < 0.) gp->delta_HI = -gp->HI_density;
	if (gp->HI_density + gp->delta_HI > gp->rhoH) gp->delta_HI = gp->rhoH - gp->HI_density;
	if (gp->HII_density + gp->delta_HII < 0.) gp->delta_HII = -gp->HII_density;
	if (gp->HII_density + gp->delta_HII > gp->rhoH) gp->delta_HII = gp->rhoH - gp->HII_density;
	if (gp->H2I_density + gp->delta_H2I < 0.) gp->delta_H2I = -gp->H2I_density;
	if (gp->H2I_density + gp->delta_H2I > gp->rhoH) gp->delta_H2I = gp->rhoH - gp->H2I_density;
	if (gp->H2II_density + gp->delta_H2II < 0.) gp->delta_H2II = -gp->H2II_density;
	if (gp->H2II_density + gp->delta_H2II > gp->rhoH) gp->delta_H2II = gp->rhoH - gp->H2II_density;
	if (gp->HM_density + gp->delta_HM < 0.) gp->delta_HM = -gp->HM_density;
	if (gp->HM_density + gp->delta_HM > gp->rhoH) gp->delta_HM = gp->rhoH - gp->HM_density;

	/* Evolve H species */
	gp->HI_density += gp->delta_HI;
	gp->HII_density += gp->delta_HII;
	gp->H2I_density += gp->delta_H2I;
	gp->H2II_density += gp->delta_H2II;
	gp->HM_density += gp->delta_HM;

	/* Evolve HeI, HeII */
	gp->HeI_density += gp->delta_HeI;
	gp->HeII_density += gp->delta_HeII;
	gp->HeIII_density += gp->delta_HeIII;
	if (gp->HeI_density < tiny) gp->HeI_density = tiny;
	if (gp->HeII_density < tiny) gp->HeII_density = tiny;
	if (gp->HeIII_density < tiny) gp->HeIII_density = tiny;

	/* Ensure H & He are conserved (since these are unaffected by dust) */
	const double H_old = gp_old->HI_density + gp_old->HII_density + gp_old->H2I_density + gp_old->H2II_density + gp_old->HM_density;
	const double He_old = gp_old->HeI_density + gp_old->HeII_density + gp_old->HeIII_density;
	const double H_new = gp->HI_density + gp->HII_density + gp->H2I_density + gp->H2II_density + gp->HM_density;
	const double He_new = gp->HeI_density + gp->HeII_density + gp->HeIII_density;

	const double H_factor = H_old / (H_new+tiny);
	const double He_factor = He_old / (He_new+tiny);

	/* Correct densities to conserve H, He */
	gp->HI_density *= H_factor;
	gp->HII_density *= H_factor;
	gp->H2I_density *= H_factor;
	gp->H2II_density *= H_factor;
	gp->HM_density *= H_factor;

	gp->HeI_density *= He_factor;
	gp->HeII_density *= He_factor;
	gp->HeIII_density *= He_factor;

	check_hydrogen(gp, 2);

	/* Compute new electron density */
	compute_electron_density(gp);

	return;
}

void evolve_pred_corr(grackle_part_data *gp, grackle_part_data *gp_old, chemistry_data *chemistry)
{
	/* Evolve all primordial species */
	gp->HI_density = 0.5 * (gp->HI_density + gp_old->HI_density);
	gp->HII_density = 0.5 * (gp->HII_density + gp_old->HII_density);
	gp->H2I_density = 0.5 * (gp->H2I_density + gp_old->H2I_density);
	gp->H2II_density = 0.5 * (gp->H2II_density + gp_old->H2II_density);
	gp->HM_density = 0.5 * (gp->HM_density + gp_old->HM_density);

	gp->HeI_density = 0.5 * (gp->HeI_density + gp_old->HeI_density);
	gp->HeII_density = 0.5 * (gp->HeII_density + gp_old->HeII_density);
	gp->HeIII_density = 0.5 * (gp->HeIII_density + gp_old->HeIII_density);

	/* Now ensure H & He are conserved */
	const double H_old = gp_old->HI_density + gp_old->HII_density + gp_old->H2I_density + gp_old->H2II_density + gp_old->HM_density;
	const double He_old = gp_old->HeI_density + gp_old->HeII_density + gp_old->HeIII_density;
	const double H_new = gp->HI_density + gp->HII_density + gp->H2I_density + gp->H2II_density + gp->HM_density;
	const double He_new = gp->HeI_density + gp->HeII_density + gp->HeIII_density;

	const double H_factor = H_old / (H_new+tiny);
	const double He_factor = He_old / (He_new+tiny);

	/* Correct densities proportionally to conserve H, He */
	gp->HI_density *= H_factor;
	gp->HII_density *= H_factor;
	gp->H2I_density *= H_factor;
	gp->H2II_density *= H_factor;
	gp->HM_density *= H_factor;

	gp->HeI_density *= He_factor;
	gp->HeII_density *= He_factor;
	gp->HeIII_density *= He_factor;

	/* Compute new electron density */
	compute_electron_density(gp);

	return;
}

void save_species_data(grackle_part_data *gp_orig, grackle_part_data *gp_target) 
{
        /* Copy species fractions from gp_orig to gp_target */
	gp_target->HI_density = gp_orig->HI_density;
	gp_target->HII_density = gp_orig->HII_density;
	gp_target->H2I_density = gp_orig->H2I_density;
	gp_target->H2II_density = gp_orig->H2II_density;
	gp_target->HM_density = gp_orig->HM_density;
	gp_target->HeI_density = gp_orig->HeI_density;
	gp_target->HeII_density = gp_orig->HeII_density;
	gp_target->HeIII_density = gp_orig->HeIII_density;
	gp_target->e_density = gp_orig->e_density;
	return;
}

void crackle_cooling_time(grackle_field_data *p, chemistry_data *chemistry, chemistry_data_storage grackle_rates, photo_rate_storage my_uvb_rates, code_units *units, gr_float *tcool)
{
	grackle_part_data gp[1];
	chemistry_rate_storage my_rates;  // interpolated rates for this field
	interp_struct interpolation;
	crackle_units cunits;

	/* Copy to grackle_part_data */
	copy_grackle_fields_to_part(p, gp, chemistry);
	/* Compute some basic properties */
	set_rhot(gp, units, chemistry);
	/* Set up various unit conversions etc */
	set_crackle_units(units, grackle_rates, chemistry->Gamma, &cunits);
	/* initialize inteprolation of chemistry rate tables for this particle */
	init_temperature_interpolation(gp, chemistry, &interpolation, &cunits, &grackle_rates, &my_uvb_rates);
	/* Set up cooling/heating rates interpolation */
	setup_temperature_interpolation(gp->logtem, chemistry, &interpolation);
	/* Get interpolated chemistry rates for this particle */
	lookup_chemistry_coeffs(chemistry->primordial_chemistry, &grackle_rates, &my_rates, &interpolation);  
	/* ISM flag: see whether we should be forming dust and evolving H2 */
	int ism_flag = (gp->isrf_habing >= 0.);
	if (gp->isrf_habing < 0.) gp->isrf_habing = 0.;
	/* Compute rate of change of thermal energy */
	compute_edot(gp, chemistry, &grackle_rates, &my_rates, &my_uvb_rates, &interpolation, units, &cunits, ism_flag);  

	*tcool = gp->internal_energy * gp->density / (gp->edot+tiny);
	return;
}

void crackle_temperature(grackle_field_data *p, chemistry_data *chemistry, code_units *units, gr_float *temp)
{
	grackle_part_data gp;

	/* Copy to grackle_part_data */
	copy_grackle_fields_to_part(p, &gp, chemistry);
	/* Compute some basic properties */
	set_rhot(&gp, units, chemistry);

	*temp = gp.tgas;
	return;
}

