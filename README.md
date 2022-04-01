# icl-hbv
The HBV model developed at Imperial College London, as described in the paper ["The impact of the timely birth dose vaccine on the global elimination of hepatitis B"](https://www.nature.com/articles/s41467-021-26475-6).

All of the scripts were run in MATLAB (version R2020b).

In the folder `src`, the folder `analysis` contains the scripts to run.

`main_script.m` runs `country_level_analyses.m`, which in turn runs `HBVmodel.m` (the latter is in the folder `model`), to generate 800 results files containing results for each country.

`summarize_stochastic_runs.m` summarizes the results from the 800 results files into results for each country, each WHO region, as well as global results.

`draw_figures.m` draws Figures 1e&ndash;1f, 2a&ndash;2d, 3a, 3b, 4a&ndash;4b, 5a&ndash;5b, Supplementary figures 1, 2a&ndash;2b and 3 and Supplementary Tables 3 and 4 (default analyses). Data structures for running this code can be downloaded from the Dropbox folder <https://www.dropbox.com/sh/veie946cwcl4y4m/AAA7X22L22cvqTC-tcm4HzJga?dl=0> and put into the folder `outputs`.

The folder `data_and_script_for_figures` contains the script `draw_figs.m` for drawing Figures 1f, 2b, 3a, 3b, 4a&ndash;4b, 5a&ndash;5b, Supplementary Figures 1, 2b and 3 and Supplementary Tables 3 and 4 (default analyses). The data structure for running this code can be downloaded from the Dropbox folder <https://www.dropbox.com/sh/veie946cwcl4y4m/AAA7X22L22cvqTC-tcm4HzJga?dl=0> and put into the folder `data_and_script_for_figures`.

