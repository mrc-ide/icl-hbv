# icl-hbv
A version of the HBV model developed at Imperial College London for investigating the impact and cost-effectiveness of peripartum antiviral prophylaxis, as described in the paper "Impact and cost-effectiveness of HBV prophylaxis in pregnancy: a modelling study".

All of the scripts were run in MATLAB (version R2022a).

In the folder `src`, the folder `analysis` contains the scripts to run.

`main_script.m` runs `country_level_analyses.m`, which in turn runs `HBVmodel.m` (the latter is in the folder `model`), to generate results files containing results for each country.

`make_cost_objects_countries.m` and `cost_model_countries.m` (the latter also in the folder `model`) contain the cost assumptions and calculations, respectively.

