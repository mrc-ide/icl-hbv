# icl-hbv

Scripts for generating country-level disease burden estimates using the HBV model developed at Imperial College London, as described in the article ["Requirements for global elimination of hepatitis B: a modelling study"](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(16)30204-3/fulltext) and updated in the article ["The impact of the timely birth dose vaccine on the global elimination of hepatitis B"](https://www.nature.com/articles/s41467-021-26475-6).

All of the scripts were run in MATLAB (version R2023a).

In the folder `src`, the folder `analysis` contains the scripts to run and the folder `model` contains model files called by the scripts in `analysis`. Some of the data structures for running the scripts are in the folder `resources` (vaccination coverage data has been withheld), and the scripts save their results to the folder `outputs`.

The script `main_script_prep_data.m` contains functions to prepare data structures containing deomographic and epidemiolgical data for use in subsequent scripts:
- The functions `make_regional_map`, `import_demographic_data_Montagu`, `process_demographic_data_Montagu` and `make_parameters_map` together generate a cell array of the ISO codes `ListOfAllISOs.mat` as well as `regional_data_GAVI_run.mat` and `params_GAVI_run.mat`, which contain country-level demographic and epidemiolgical data as well as disability weights.
- The function `import_coverage_data_Montagu` generates `scenario_lists_GAVI_run.mat`, which contains country-level coverage data for the various vaccination scenarios.
- The functions `make_cda_ott_HBsAg_HBeAg_prevalence_maps` and `make_articles_data_map` generate `reference_prevalences.mat` and `articles_data_map.mat`, which contain the country-level HBsAg and HBeAg/HBsAg prevelences used for model calibration and initialisation.
- The function `make_CDA_treatment_2016_map` generates `CDA_treatment_2016.mat`, which contains the number of individiuals in treatment in countries in 2016.
- The functions `make_GBD_HCCdeaths_data_map`, `make_GBD_HCC_pop_sizes_data_map`, `make_GBD_cirrh_deaths_data_map` and `make_GBD_cirrh_pop_sizes_data_map` generate `GBD_cirrh_deaths_1990_Number_data_202310_map.mat`, `GBD_cirrh_deaths_1990_Rate_data_202310_map.mat`, `GBD_cirrh_deaths_2005_Number_data_202310_map.mat`, `GBD_cirrh_deaths_2005_Rate_data_202310_map.mat`, `GBD_cirrh_popSizes_1990_data_202310_map.mat`, `GBD_cirrh_popSizes_2005_data_202310_map.mat`, `GBD_HCCdeaths_1990_Number_data_202310_map.mat`, `GBD_HCCdeaths_1990_Rate_data_202310_map.mat`, `GBD_HCCdeaths_2005_Number_data_202310_map.mat`, `GBD_HCCdeaths_2005_Rate_data_202310_map.mat`, `GBD_HCCpopSizes_1990_data_202310_map.mat` and `GBD_HCCpopSizes_2005_data_202310_map.mat`, which contain country-level HBV-related deaths data.

The script `main_script_optimize_7_params.m`, which runs `HBVmodel_optimize.m` (the latter is in the folder `model`), calibrates the model to country-level HBV prevalence and HBV-related deaths data using the ABC SMC algorithm.
The script `fitting_plots_code.m` draws plots of the fit for each country.
The script `prepare_data_structures_from_calibrations.m` generates `country_optimization_s_e_HCCdeaths.mat`, which contains the country-level HBsAg and HBeAg/HBsAg prevelences used for model calibration and initialisation, and `country_theta_hats.mat`, which contains the 200 particle values for running the model for each country.
The script `make_stocas_params_mat.m` generates `params_GAVI_run_stochastic.mat`, which contains the 200 particle values for running the model for each country.

The script `make_pop_size_HBsAg_treatment_map.m` generates `pop_size_HBsAg_treatment_2016_map.mat`, which contains country-level treatment data in 2016.
The script `main_script_determine_treatment_rates_countries.m`, which runs `HBVmodel_treatment.m` (the latter is in the folder `model`), finds the treatment rate in a country that results in the desired level of treatment amongst the treatment-eligible by 2030 (either the same level of treatment in 2030 as in 2016 or 40% of treatment-eligible in treatment by 2030).
The script `make_treatment_rates_mat_from_stochastic_runs.m` generates `treatment_lower_upper_rates_particles_1.mat` and `treatment_lower_upper_rates_particles_25.mat`, which contain the treatment rates to be used in the model runs for each country.

The scripts `measure_mean_durations_particles_200.m` and `make_mean_times_mat_from_stochastic_runs.m` generate `mean_median_duration_times_mort_migr_treat_particles_200.mat`, which contains country-level estimates of the mean length of time spent in the disease states for a cohort of a particular sex and age starting from the year 2000.

The script `main_script_mont_stochas.m`, which runs `HBVmodel_Montagu.m` (the latter is in the folder `model`), generates country-level estimates of the model outcomes.

The scripts `make_geometric_median_global.m` and `Montagu_stochas_comparisons.m` add up the country-level estimates to generate global-level estimates of the model outcomes.

The scripts `compare_geo_med_model_outputs_of_scenarios.m` and `compare_summarized_model_outputs_of_scenarios.m` check that the country-level results from the various vaccination scenarios make sense relative to one another.

The scripts `write_to_csv_geo_med.m`, `write_to_csv_means.m` and `write_to_csv_stochastic.m` generate csv files of the results of the model runs resolved by country, year, age and model outcome.


