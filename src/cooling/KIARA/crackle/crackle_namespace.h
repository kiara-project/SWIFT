/*******************************************************************************
 * This file is part of SWIFT.
 *
 * crackle namespace: every symbol exported by the crackle library is renamed
 * with a crk_ prefix (names already starting with crackle_ are kept), so that
 * it cannot clash with, or be interposed by, an external Grackle library.
 * This file must be included before any other crackle/Grackle header; the
 * crackle sources are compiled with -include crackle_namespace.h and SWIFT
 * gets it through crackle/grackle.h.
 *
 * Generated from the exported symbols of the imported library (nm).
 ******************************************************************************/
#ifndef CRACKLE_NAMESPACE_H
#define CRACKLE_NAMESPACE_H

#define add_h2dust_reaction_rate crk_add_h2dust_reaction_rate
#define GAH2_rate crk_GAH2_rate
#define GAHI_rate crk_GAHI_rate
#define GAHe_rate crk_GAHe_rate
#define GAHp_rate crk_GAHp_rate
#define GAel_rate crk_GAel_rate
#define GP99HighDensityLimit_rate crk_GP99HighDensityLimit_rate
#define GP99LowDensityLimit_rate crk_GP99LowDensityLimit_rate
#define H2LTE_rate crk_H2LTE_rate
#define HDlow_rate crk_HDlow_rate
#define HDlte_rate crk_HDlte_rate
#define _build_index_helper crk__build_index_helper
#define _calculate_cooling_time crk__calculate_cooling_time
#define _calculate_gamma crk__calculate_gamma
#define _calculate_pressure crk__calculate_pressure
#define _calculate_temperature crk__calculate_temperature
#define _free_cloudy_data crk__free_cloudy_data
#define _k13dd_rate crk__k13dd_rate
#define _solve_chemistry crk__solve_chemistry
#define add_h2dust_reaction_rate crk_add_h2dust_reaction_rate
#define add_k13dd_reaction_rate crk_add_k13dd_reaction_rate
#define add_reaction_rate crk_add_reaction_rate
#define add_scalar_reaction_rate crk_add_scalar_reaction_rate
#define brem_rate crk_brem_rate
#define calculate_cooling_time crk_calculate_cooling_time
#define calculate_gamma crk_calculate_gamma
#define calculate_pressure crk_calculate_pressure
#define calculate_temperature crk_calculate_temperature
#define ceHI_rate crk_ceHI_rate
#define ceHeII_rate crk_ceHeII_rate
#define ceHeI_rate crk_ceHeI_rate
#define ciHI_rate crk_ciHI_rate
#define ciHeII_rate crk_ciHeII_rate
#define ciHeIS_rate crk_ciHeIS_rate
#define ciHeI_rate crk_ciHeI_rate
#define cie_table_c crk_cie_table_c
#define cie_thin_cooling_rate crk_cie_thin_cooling_rate
#define cieco_rate crk_cieco_rate
#define comp_rate crk_comp_rate
#define compute_HIdot crk_compute_HIdot
#define compute_dedot crk_compute_dedot
#define compute_edot crk_compute_edot
#define compute_iteration_dt crk_compute_iteration_dt
#define evolve_H2 crk_evolve_H2
#define evolve_elements crk_evolve_elements
#define evolve_energy_ros2 crk_evolve_energy_ros2
#define evolve_helium crk_evolve_helium
#define evolve_hydrogen crk_evolve_hydrogen
#define evolve_internal_energy crk_evolve_internal_energy
#define evolve_pred_corr crk_evolve_pred_corr
#define free_chemistry_data crk_free_chemistry_data
#define gamma_isrf_rate crk_gamma_isrf_rate
#define gammah_rate crk_gammah_rate
#define gasGrain_rate crk_gasGrain_rate
#define get_temperature_units crk_get_temperature_units
#define get_velocity_units crk_get_velocity_units
#define grackle_data crk_grackle_data
#define grackle_num_params crk_grackle_num_params
#define grackle_rates crk_grackle_rates
#define grackle_verbose crk_grackle_verbose
#define h2dust_rate crk_h2dust_rate
#define h2k01_rate crk_h2k01_rate
#define hyd01k_rate crk_hyd01k_rate
#define initialize_UVbackground_data crk_initialize_UVbackground_data
#define initialize_chemistry_data crk_initialize_chemistry_data
#define initialize_cloudy_data crk_initialize_cloudy_data
#define initialize_empty_UVBtable_struct crk_initialize_empty_UVBtable_struct
#define initialize_empty_chemistry_data_storage_struct crk_initialize_empty_chemistry_data_storage_struct
#define initialize_empty_cloudy_data_struct crk_initialize_empty_cloudy_data_struct
#define initialize_rates crk_initialize_rates
#define k10_rate crk_k10_rate
#define k11_rate crk_k11_rate
#define k12_rate crk_k12_rate
#define k13_rate crk_k13_rate
#define k13dd_rate crk_k13dd_rate
#define k14_rate crk_k14_rate
#define k15_rate crk_k15_rate
#define k16_rate crk_k16_rate
#define k17_rate crk_k17_rate
#define k18_rate crk_k18_rate
#define k19_rate crk_k19_rate
#define k1_rate crk_k1_rate
#define k20_rate crk_k20_rate
#define k21_rate crk_k21_rate
#define k22_rate crk_k22_rate
#define k23_rate crk_k23_rate
#define k2_rate crk_k2_rate
#define k3_rate crk_k3_rate
#define k4_rate crk_k4_rate
#define k50_rate crk_k50_rate
#define k51_rate crk_k51_rate
#define k52_rate crk_k52_rate
#define k53_rate crk_k53_rate
#define k54_rate crk_k54_rate
#define k55_rate crk_k55_rate
#define k56_rate crk_k56_rate
#define k57_rate crk_k57_rate
#define k58_rate crk_k58_rate
#define k5_rate crk_k5_rate
#define k6_rate crk_k6_rate
#define k7_rate crk_k7_rate
#define k8_rate crk_k8_rate
#define k9_rate crk_k9_rate
#define local_calculate_cooling_time crk_local_calculate_cooling_time
#define local_calculate_gamma crk_local_calculate_gamma
#define local_calculate_pressure crk_local_calculate_pressure
#define local_calculate_temperature crk_local_calculate_temperature
#define local_calculate_temperature_table crk_local_calculate_temperature_table
#define local_chemistry_data_access_double crk_local_chemistry_data_access_double
#define local_chemistry_data_access_int crk_local_chemistry_data_access_int
#define local_chemistry_data_access_string crk_local_chemistry_data_access_string
#define local_free_chemistry_data crk_local_free_chemistry_data
#define local_initialize_chemistry_data crk_local_initialize_chemistry_data
#define local_initialize_chemistry_parameters crk_local_initialize_chemistry_parameters
#define local_solve_chemistry crk_local_solve_chemistry
#define logT_spacing crk_logT_spacing
#define logT_spacing_dust crk_logT_spacing_dust
#define n_cr_d1_rate crk_n_cr_d1_rate
#define n_cr_d2_rate crk_n_cr_d2_rate
#define n_cr_n_rate crk_n_cr_n_rate
#define param_name_double crk_param_name_double
#define param_name_int crk_param_name_int
#define param_name_string crk_param_name_string
#define reHII_rate crk_reHII_rate
#define reHeII1_rate crk_reHeII1_rate
#define reHeII2_rate crk_reHeII2_rate
#define reHeIII_rate crk_reHeIII_rate
#define read_dataset crk_read_dataset
#define regr_rate crk_regr_rate
#define roth_rate crk_roth_rate
#define rotl_rate crk_rotl_rate
#define save_species_data crk_save_species_data
#define set_default_chemistry_parameters crk_set_default_chemistry_parameters
#define set_velocity_units crk_set_velocity_units
#define show_parameters crk_show_parameters
#define solve_chemistry crk_solve_chemistry
#define t_cie_c crk_t_cie_c
#define update_UVbackground_rates crk_update_UVbackground_rates
#define vibh_rate crk_vibh_rate

#endif /* CRACKLE_NAMESPACE_H */
