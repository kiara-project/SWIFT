
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DUST_FRAC_CAP 0.9 // Total dust growth is limited to this fraction of gas-phase metals
#define DUST_CHANGE_HI 100.0  // factor by which dust density can go up in a single step
#define DUST_CHANGE_LO 0.5  // factor by which dust density can go down in a single step
#define T_SPUTTERING 2.e4  // Temeprature below which sputtering is negligible
#define DUST_GROWTH_TAU_GYR 1.0  // dust growth time-scale (Gyr) at density dust_growth_densref

static inline void evolve_dust(grackle_part_data *gp, chemistry_data *chemistry, code_units *units, int ism_flag, double dtit)
{
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
                //const double Ms100 = 6800.0*chemistry->sne_coeff*(100.0/chemistry->sne_shockspeed)*(100.0/chemistry->sne_shockspeed) * SolarMass / mass_unit;  //code units, gas mass shocked per SNe (Sedov-Taylor phase) in Simba
		const double SN_dust_destr_rate = 22.;  // dust destruction rate per SN in Mo/Myr (Kirchschlager+24, model BH)
		const double tdyn = 1./sqrt(6.67e-8 * dens_cgs) / (1.e6 * sec_per_year);  // assume SN remnant lasts for tdyn; convert to Myr
		const double Ms100 = SN_dust_destr_rate * tdyn * SolarMass / (mass_unit * chemistry->dust_destruction_eff);  // dust mass destroyed per SN, in code units
                f_shocked = min((Ms100 * gp->SNe_density) / rhogas, 1.);  // fraction of mass shock-heated
                drhos = gp->dust_density * f_shocked * chemistry->dust_destruction_eff;  // some fraction of dust is destroyed in shocked gas
	//if (gp->verbose) printf("dust: %g %g %g %g %g %g %g\n",gp->dust_density, Ms100 * mass_unit * chemistry->dust_destruction_eff / SolarMass, tdyn, dens_cgs, f_shocked, drhos, (gp->metal_density/(gp->metal_density+gp->dust_density)) * (gp->dust_density/tau_accr) * dtit);
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

