# icl-hbv
The HBV model developed at Imperial College London.

All of the scripts run in MATLAB (version R2020b).

main_script.m runs country_level_analyses.m, which in turn runs HBVmodel.m, to generate 800 results files containing results for each country.

summarize_stochastic_runs.m summarizes the results from the 800 files into results for each country, each WHO region, as well as global results.

The folder data_and_script_for_figures contains the data (regions_countries_data.mat) and script (draw_figs.m) for drawing figures 1f, 2b, 3a, 3b, 4a, 4b, 5a and 5b, Supplementary figures 1, 2b and 3 and Supplementary Tables 3 and 4 (default analyses).

