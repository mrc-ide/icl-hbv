%%

% This script contains code for drawing 
% Figures 1e--1f, 2a--2d, 3a, 3b, 4a--4b, 5a--5b, 
% Supplementary Figures 1 and 2a--2b and
% Supplementary Tables 3 and 4 (default analyses).


clear


sensitivity_analysis_list = {'default','infant_100','treat_medium','treat_high'};
num_sensitivity_analyses = length(sensitivity_analysis_list);


num_stochas_runs = 200;


Regions_map = load('WHO_region_map.mat').WHO_region_map;
ListOfScenarioLabels = load('scenarios_array.mat').label_array;
filename = 'results_countries_default_stochastic_run_1.mat';
results_default_run_1_map = load(filename).outMap;


num_scenarios = length(ListOfScenarioLabels);
assert(num_scenarios==12)
ListOfScenarios = {...
    'Status quo infant & BD',...
    'Status quo infant & BD expansion to 25%',...
    'Status quo infant & BD expansion to 50%',...
    'Status quo infant & BD expansion to 75%',...
    'Status quo infant & BD expansion to 90%',...
    'Status quo infant & BD drop 5 2020',...
    'Status quo infant & BD drop 10 2020',...
    'Status quo infant & BD drop 15 2020',...
    'Status quo infant & BD drop 20 2020',...
    'Status quo infant & BD delayed expansion 2023 to 2030',...
    'Status quo infant & BD delayed expansion 2023 to 2033',...
    'Status quo infant & BD delayed expansion 2025 to 2040'};
assert(length(ListOfScenarios)==num_scenarios)

scenario_nums_struct.status_quo_infant_status_quo_BD_num = 1;
scenario_nums_struct.status_quo_infant_BD_to_25_num = 2;
scenario_nums_struct.status_quo_infant_BD_to_50_num = 3;
scenario_nums_struct.status_quo_infant_BD_to_75_num = 4;
scenario_nums_struct.status_quo_infant_BD_to_90_num = 5;
scenario_nums_struct.status_quo_infant_BD_drop_5_num = 6;
scenario_nums_struct.status_quo_infant_BD_drop_10_num = 7;
scenario_nums_struct.status_quo_infant_BD_drop_15_num = 8;
scenario_nums_struct.status_quo_infant_BD_drop_20_num = 9;
scenario_nums_struct.BD_delayed_2023_2030_num = 10;
scenario_nums_struct.BD_delayed_slowed_2023_2033_num = 11;
scenario_nums_struct.BD_delayed_slowed_2025_2040_num = 12;
assert(length(fieldnames(scenario_nums_struct))==num_scenarios)


infilepath_stochastic_results = '';
file_string = 'default';
[regions_cell_array_default,countries_cell_array_default] = make_stochas_regions_countries_cell_arrays(...
    infilepath_stochastic_results,file_string,...
    num_scenarios,scenario_nums_struct);
disp(['Finished making data structures for the default sensitivity analysis.'])


results_stochas_regions_StatusQuo_map = regions_cell_array_default{scenario_nums_struct.status_quo_infant_status_quo_BD_num};
results_stochas_countries_StatusQuo_map = countries_cell_array_default{scenario_nums_struct.status_quo_infant_status_quo_BD_num};


regions_global_list = keys(results_stochas_regions_StatusQuo_map);
num_regions_global = length(regions_global_list);
assert(num_regions_global==6+1)
regions_list = regions_global_list;
index_global = strcmp(regions_list,'Global');
regions_list(index_global) = '';
num_regions = length(regions_list);
assert(num_regions==6)
index_global = find(index_global); % the index number


countries_list = keys(results_stochas_countries_StatusQuo_map);
num_countries = length(countries_list);
assert(num_countries==110)


results_stochas_Global_StatusQuo_struct = results_stochas_regions_StatusQuo_map('Global');

calandar_years_vec = 1980:2100;
num_calandar_years = length(calandar_years_vec);
ilower = find(results_stochas_Global_StatusQuo_struct.Time>=calandar_years_vec(1), 1);
iupper = find(results_stochas_Global_StatusQuo_struct.Time>=calandar_years_vec(end), 1);
calandar_years_pos_vec = ilower:iupper;
i2025 = find(results_stochas_Global_StatusQuo_struct.Time>=2025, 1);

cohort_years_vec = 2015:2050;
num_cohort_years = length(cohort_years_vec);
ilower = find(results_stochas_Global_StatusQuo_struct.Time>=cohort_years_vec(1), 1);
iupper = find(results_stochas_Global_StatusQuo_struct.Time>=cohort_years_vec(end), 1);
cohort_years_pos_vec = ilower:iupper;

num_year_divisions = 10;
dt = 1/num_year_divisions;
num_time_steps = length(1980:dt:2100);


vaccination_map_default = containers.Map;
for scenario_num=1:num_scenarios

    scenario = ListOfScenarios{scenario_num};
    countryMap = results_default_run_1_map(scenario);

    scenario_map_default = containers.Map;
    for country_num=1:num_countries

        ISO = countries_list{country_num};
        scenario_run_1_country_struct = countryMap(ISO);
        scenario_country_struct.InfantVacc = scenario_run_1_country_struct.InfantVacc;
        scenario_country_struct.BirthDoseVacc = scenario_run_1_country_struct.BirthDoseVacc;
        scenario_country_struct.Tot_Pop_2025 = scenario_run_1_country_struct.Tot_Pop_1yr(:,i2025);

        scenario_map_default(ISO) = scenario_country_struct;

    end

    vaccination_map_default(scenario) = scenario_map_default;

end % end scenario_num for loop


region_colour_mat = [...
    0 0 255; ...     % blue    % AFRO
    0 204 0; ...     % green   % EMRO
    255 128 0; ...   % orange  % EURO
    128 128 128; ... % grey    % PAHO
    255 51 255; ...  % magenta % SEARO
    153 76 0 ...     % brown   % WPRO
    ]/255;
assert(isequal(size(region_colour_mat),[num_regions 3]))


scenario_line_style_array = {'-','-',':','-','-','-','-','-','--','-',':','--'};
scenario_colour_mat = [...
    0 0 0; ...       % black   % 1
    204 0 0; ...     % ???     % 2a
    255 0 0; ...     % red     % 2b
    255 102 102; ... % ???     % 2c
    255 153 153; ... % ???     % 2d
    204 102 0; ...   % brown   % 3a
    255 153 51; ...  % orange  % 3b
    255 255 0; ...   % yellow  % 3c
    192 192 192; ... % grey    % 3d
    0 153 0; ...     % ???     % 4a
    0 204 0; ...     % ???     % 4b
    0 255 0 ...      % green   % 4c
    ]/255;
assert(length(scenario_line_style_array)==num_scenarios)
assert(isequal(size(scenario_colour_mat),[num_scenarios 3]))






country_region_cell_array = draw_Suppl_Fig_1_fun(Regions_map,countries_list);
% The cell array country_region_cell_array can be used to draw a map in Excel for Supplementary Fig. 1.

draw_Figs_1e_1f_Suppl_Figs_2a_2b_fun(...
    num_year_divisions,...
    calandar_years_vec,num_calandar_years,...
    ListOfScenarios,ListOfScenarioLabels,scenario_nums_struct,scenario_line_style_array,scenario_colour_mat,...
    num_regions,regions_list,num_regions_global,Regions_map,region_colour_mat,...
    num_countries,countries_list,...
    vaccination_map_default)

draw_Figs_2a_2b_2c_2d_fun(...
    num_stochas_runs,...
    calandar_years_vec,...
    ListOfScenarios,ListOfScenarioLabels,scenario_nums_struct,scenario_line_style_array,scenario_colour_mat,...
    countries_list,...
    regions_cell_array_default)

draw_Fig_3a_fun(...
    num_stochas_runs,...
    cohort_years_vec,num_cohort_years,cohort_years_pos_vec,...
    scenario_nums_struct,...
    num_regions,regions_list,region_colour_mat,...
    countries_list,...
    regions_cell_array_default)

[countries_list, perc_global_total_excess_deaths_vec, perc_global_reduc_in_deaths_vec] = draw_Fig_3b_fun(...
    num_stochas_runs,...
    cohort_years_vec,cohort_years_pos_vec,...
    scenario_nums_struct,...
    num_countries,countries_list,...
    countries_cell_array_default);
% This code generates the vector perc_global_total_excess_deaths_vec that, together with the cell array countries_list, can be used to draw a map in Excel for Fig. 3b.
% This code also generates the vector perc_global_reduc_in_deaths_vec that, together with the cell array countries_list, can be used to draw a map in Excel for Supplementary Fig. 3.

Fig_4a_data_struc = prepare_Fig_4a_data_struc_fun(...
    num_stochas_runs,...
    calandar_years_vec,calandar_years_pos_vec,...
    ListOfScenarioLabels,scenario_nums_struct,scenario_colour_mat,...
    num_regions,regions_list,num_regions_global,regions_global_list,index_global,...
    countries_list,...
    regions_cell_array_default);

Fig_4b_data_struc = prepare_Fig_4b_data_struc_fun(...
    num_stochas_runs,num_year_divisions,...
    calandar_years_vec,num_calandar_years,calandar_years_pos_vec,...
    ListOfScenarios,scenario_nums_struct,...
    num_regions,regions_list,Regions_map,region_colour_mat,...
    num_countries,countries_list,...
    countries_cell_array_default,vaccination_map_default);

draw_Figs_4a_4b_fun(Fig_4a_data_struc,Fig_4b_data_struc)

Fig_5a_data_struc = prepare_Fig_5a_data_struc_fun(...
    num_stochas_runs,...
    cohort_years_vec,cohort_years_pos_vec,...
    scenario_nums_struct,...
    num_regions,regions_list,num_regions_global,regions_global_list,index_global,region_colour_mat,... 
    countries_list,...
    regions_cell_array_default);

Fig_5b_data_struc = prepare_Fig_5b_data_struc_fun(...
    num_stochas_runs,...
    cohort_years_vec,cohort_years_pos_vec,...
    scenario_nums_struct,...
    num_regions,regions_list,num_regions_global,regions_global_list,index_global,region_colour_mat,...
    countries_list,...
    regions_cell_array_default);

draw_Figs_5a_5b_fun(Fig_5a_data_struc,Fig_5b_data_struc)

draw_Suppl_Tables_3_4_fun(...
    num_stochas_runs,...
    calandar_years_vec,calandar_years_pos_vec,...
    cohort_years_vec,cohort_years_pos_vec,...
    scenario_nums_struct,...
    num_regions_global,regions_global_list,...
    countries_list,...
    regions_cell_array_default)





function country_region_cell_array = draw_Suppl_Fig_1_fun(Regions_map,countries_list)

    country_region_cell_array = countries_list';
    country_region_cell_array(:,2) = cellfun(@(xx) Regions_map(xx),country_region_cell_array(:,1),'UniformOutput',false);
    disp("Use the cell array country_region_cell_array to draw a map in Excel.")

end % end of draw_Suppl_Fig_1_fun function





function draw_Figs_1e_1f_Suppl_Figs_2a_2b_fun(...
    num_year_divisions,...
    calandar_years_vec,num_calandar_years,...
    ListOfScenarios,ListOfScenarioLabels,scenario_nums_struct,scenario_line_style_array,scenario_colour_mat,...
    num_regions,regions_list,num_regions_global,Regions_map,region_colour_mat,...
    num_countries,countries_list,...
    vaccination_map_default)


    years_vec_01yr = calandar_years_vec(1):(1/num_year_divisions):calandar_years_vec(end);
    num_time_steps = length(years_vec_01yr);


    status_quo_infant_status_quo_BD_num = scenario_nums_struct.status_quo_infant_status_quo_BD_num;
    status_quo_infant_BD_to_90_num = scenario_nums_struct.status_quo_infant_BD_to_90_num;
    status_quo_infant_BD_drop_20_num = scenario_nums_struct.status_quo_infant_BD_drop_20_num;
    BD_delayed_slowed_2025_2040_num = scenario_nums_struct.BD_delayed_slowed_2025_2040_num;


    [countries_infant_coverage_1980_2100_status_quo_mat, ...
        countries_infant_coverage_1980_2100_BD_90_mat, ...
        countries_infant_coverage_1980_2100_BD_drop_20_mat, ...
        countries_infant_coverage_1980_2100_BD_2025_2040_mat] = ...
        deal(zeros(num_countries,num_time_steps));
    [countries_BD_coverage_1980_2100_status_quo_mat, ...
        countries_BD_coverage_1980_2100_BD_90_mat, ...
        countries_BD_coverage_1980_2100_BD_drop_20_mat, ...
        countries_BD_coverage_1980_2100_BD_2025_2040_mat] = ...
        deal(zeros(num_countries,num_time_steps));
    [countries_Tot_Pop_2025_status_quo_mat, ...
        countries_Tot_Pop_2025_BD_90_mat, ...
        countries_Tot_Pop_2025_BD_drop_20_mat, ...
        countries_Tot_Pop_2025_BD_2025_2040_mat] = ...
        deal(zeros(num_countries,1));


    for country_num=1:num_countries

    
        ISO = countries_list{country_num};

        [countries_infant_coverage_1980_2100_status_quo_mat,countries_BD_coverage_1980_2100_status_quo_mat,countries_Tot_Pop_2025_status_quo_mat] = ...
            get_infant_BD_Tot_Pop_2025(num_countries,...
            num_time_steps,...
            'countries',ISO,country_num,vaccination_map_default,...
            ListOfScenarios,'Status quo infant & BD',status_quo_infant_status_quo_BD_num,...
            countries_infant_coverage_1980_2100_status_quo_mat,countries_BD_coverage_1980_2100_status_quo_mat,countries_Tot_Pop_2025_status_quo_mat);

        [countries_infant_coverage_1980_2100_BD_90_mat,countries_BD_coverage_1980_2100_BD_90_mat,countries_Tot_Pop_2025_BD_90_mat] = ...
            get_infant_BD_Tot_Pop_2025(num_countries,...
            num_time_steps,...
            'countries',ISO,country_num,vaccination_map_default,...
            ListOfScenarios,'Status quo infant & BD expansion to 90%',status_quo_infant_BD_to_90_num,...
            countries_infant_coverage_1980_2100_BD_90_mat,countries_BD_coverage_1980_2100_BD_90_mat,countries_Tot_Pop_2025_BD_90_mat);

        [countries_infant_coverage_1980_2100_BD_drop_20_mat,countries_BD_coverage_1980_2100_BD_drop_20_mat,countries_Tot_Pop_2025_BD_drop_20_mat] = ...
            get_infant_BD_Tot_Pop_2025(num_countries,...
            num_time_steps,...
            'countries',ISO,country_num,vaccination_map_default,...
            ListOfScenarios,'Status quo infant & BD drop 20 2020',status_quo_infant_BD_drop_20_num,...
            countries_infant_coverage_1980_2100_BD_drop_20_mat,countries_BD_coverage_1980_2100_BD_drop_20_mat,countries_Tot_Pop_2025_BD_drop_20_mat);

        [countries_infant_coverage_1980_2100_BD_2025_2040_mat,countries_BD_coverage_1980_2100_BD_2025_2040_mat,countries_Tot_Pop_2025_BD_2025_2040_mat] = ...
            get_infant_BD_Tot_Pop_2025(num_countries,...
            num_time_steps,...
            'countries',ISO,country_num,vaccination_map_default,...
            ListOfScenarios,'Status quo infant & BD delayed expansion 2025 to 2040',BD_delayed_slowed_2025_2040_num,...
            countries_infant_coverage_1980_2100_BD_2025_2040_mat,countries_BD_coverage_1980_2100_BD_2025_2040_mat,countries_Tot_Pop_2025_BD_2025_2040_mat);

    end % end of country_num for loop


    assert(isequal(size(countries_infant_coverage_1980_2100_status_quo_mat),[num_countries num_time_steps]))
    assert(isequal(size(countries_infant_coverage_1980_2100_BD_90_mat),[num_countries num_time_steps]))
    assert(isequal(size(countries_infant_coverage_1980_2100_BD_drop_20_mat),[num_countries num_time_steps]))
    assert(isequal(size(countries_infant_coverage_1980_2100_BD_2025_2040_mat),[num_countries num_time_steps]))
    assert(all(all(countries_infant_coverage_1980_2100_status_quo_mat>=0)))
    assert(all(all(countries_infant_coverage_1980_2100_BD_90_mat>=0)))
    assert(all(all(countries_infant_coverage_1980_2100_BD_drop_20_mat>=0)))
    assert(all(all(countries_infant_coverage_1980_2100_BD_2025_2040_mat>=0)))
    assert(all(all(countries_infant_coverage_1980_2100_status_quo_mat<=100)))
    assert(all(all(countries_infant_coverage_1980_2100_BD_90_mat<=100)))
    assert(all(all(countries_infant_coverage_1980_2100_BD_drop_20_mat<=100)))
    assert(all(all(countries_infant_coverage_1980_2100_BD_2025_2040_mat<=100)))
    assert(all(all(countries_infant_coverage_1980_2100_status_quo_mat<=countries_infant_coverage_1980_2100_BD_90_mat)))
    assert(all(all(countries_infant_coverage_1980_2100_BD_drop_20_mat<=countries_infant_coverage_1980_2100_status_quo_mat)))
    assert(all(all(countries_infant_coverage_1980_2100_BD_2025_2040_mat<=countries_infant_coverage_1980_2100_BD_90_mat)))
    assert(isequal(size(countries_BD_coverage_1980_2100_status_quo_mat),[num_countries num_time_steps]))
    assert(isequal(size(countries_BD_coverage_1980_2100_BD_90_mat),[num_countries num_time_steps]))
    assert(isequal(size(countries_BD_coverage_1980_2100_BD_drop_20_mat),[num_countries num_time_steps]))
    assert(isequal(size(countries_BD_coverage_1980_2100_BD_2025_2040_mat),[num_countries num_time_steps]))
    assert(all(all(countries_BD_coverage_1980_2100_status_quo_mat>=0)))
    assert(all(all(countries_BD_coverage_1980_2100_BD_90_mat>=0)))
    assert(all(all(countries_BD_coverage_1980_2100_BD_drop_20_mat>=0)))
    assert(all(all(countries_BD_coverage_1980_2100_BD_2025_2040_mat>=0)))
    assert(all(all(countries_BD_coverage_1980_2100_status_quo_mat<=100)))
    assert(all(all(countries_BD_coverage_1980_2100_BD_90_mat<=100)))
    assert(all(all(countries_BD_coverage_1980_2100_BD_drop_20_mat<=100)))
    assert(all(all(countries_BD_coverage_1980_2100_BD_2025_2040_mat<=100)))
    assert(all(all(countries_BD_coverage_1980_2100_status_quo_mat<=countries_BD_coverage_1980_2100_BD_90_mat)))
    assert(all(all(countries_BD_coverage_1980_2100_BD_drop_20_mat<=countries_BD_coverage_1980_2100_status_quo_mat)))
    assert(all(all(countries_BD_coverage_1980_2100_BD_2025_2040_mat<=countries_BD_coverage_1980_2100_BD_90_mat)))
    assert(isequal(size(countries_Tot_Pop_2025_status_quo_mat),[num_countries 1]))
    assert(isequal(size(countries_Tot_Pop_2025_BD_90_mat),[num_countries 1]))
    assert(isequal(size(countries_Tot_Pop_2025_BD_drop_20_mat),[num_countries 1]))
    assert(isequal(size(countries_Tot_Pop_2025_BD_2025_2040_mat),[num_countries 1]))
    assert(all(all(countries_Tot_Pop_2025_status_quo_mat>=0)))
    assert(all(all(countries_Tot_Pop_2025_BD_90_mat>=0)))
    assert(all(all(countries_Tot_Pop_2025_BD_drop_20_mat>=0)))
    assert(all(all(countries_Tot_Pop_2025_BD_2025_2040_mat>=0)))


    countries_infant_coverage_1980_2100_mat = cat(3,...
        countries_infant_coverage_1980_2100_status_quo_mat,...
        countries_infant_coverage_1980_2100_BD_90_mat,...
        countries_infant_coverage_1980_2100_BD_drop_20_mat,...
        countries_infant_coverage_1980_2100_BD_2025_2040_mat); 
    countries_BD_coverage_1980_2100_mat = cat(3,...
        countries_BD_coverage_1980_2100_status_quo_mat,...
        countries_BD_coverage_1980_2100_BD_90_mat,...
        countries_BD_coverage_1980_2100_BD_drop_20_mat,...
        countries_BD_coverage_1980_2100_BD_2025_2040_mat); 
    countries_Tot_Pop_2025_mat = [...
        countries_Tot_Pop_2025_status_quo_mat,...
        countries_Tot_Pop_2025_BD_90_mat,...
        countries_Tot_Pop_2025_BD_drop_20_mat,...
        countries_Tot_Pop_2025_BD_2025_2040_mat];
    assert(isequal(size(countries_infant_coverage_1980_2100_mat),[num_countries num_time_steps 4]))
    assert(all(all(all(countries_infant_coverage_1980_2100_mat>=0))))
    assert(all(all(all(countries_infant_coverage_1980_2100_mat<=100))))
    assert(isequal(size(countries_BD_coverage_1980_2100_mat),[num_countries num_time_steps 4]))
    assert(all(all(all(countries_BD_coverage_1980_2100_mat>=0))))
    assert(all(all(all(countries_BD_coverage_1980_2100_mat<=100))))
    assert(isequal(size(countries_Tot_Pop_2025_mat),[num_countries 4]))
    assert(all(all(all(countries_Tot_Pop_2025_mat>=0))))
    countries_Tot_Pop_2025_mat = repmat(countries_Tot_Pop_2025_mat,1,1,num_time_steps);
    countries_Tot_Pop_2025_mat = permute(countries_Tot_Pop_2025_mat,[1 3 2]);
    assert(isequal(size(countries_Tot_Pop_2025_mat),[num_countries num_time_steps 4]))
    countries_infant_coverage_1980_2100_mat = countries_infant_coverage_1980_2100_mat .* countries_Tot_Pop_2025_mat; % weight timely HepB-BD coverage by population size in 2025
    assert(isequal(size(countries_infant_coverage_1980_2100_mat),[num_countries num_time_steps 4]))
    countries_BD_coverage_1980_2100_mat = countries_BD_coverage_1980_2100_mat .* countries_Tot_Pop_2025_mat; % weight timely HepB-BD coverage by population size in 2025
    assert(isequal(size(countries_BD_coverage_1980_2100_mat),[num_countries num_time_steps 4]))

    assert(isequal(size(countries_list),[1 num_countries]))
    countries_list_regions = cellfun(@(x) Regions_map(x),countries_list,'UniformOutput',false); % cell array of region of each country
    assert(isequal(size(countries_list_regions),[1 num_countries]))
    region_index_map = containers.Map(regions_list,{1,2,3,4,5,6});
    countries_index_vec = cellfun(@(x) region_index_map(x),countries_list_regions); % vector of indices from 1 to 6 for assigning a region's value to countries in that region
    assert(isequal(size(countries_index_vec),[1 num_countries]))

    regions_infant_coverage_1980_2100_mat = -99 * ones(num_regions_global,num_time_steps,4);
    regions_BD_coverage_1980_2100_mat = -99 * ones(num_regions_global,num_time_steps,4);
    regions_Tot_Pop_2025_mat = -99 * ones(num_regions_global,num_time_steps,4);
    for ii=1:num_regions_global
        if ii<=6
            region_index_vec = countries_index_vec == ii; % vector of 0s and 1s for a particular region
        else
            assert(ii==7)
            region_index_vec = ones(1,num_countries);
        end
        assert(isequal(size(region_index_vec),[1 num_countries]))

        region_index_mat = repmat(region_index_vec',1,num_time_steps,4);
        assert(isequal(size(region_index_mat),[num_countries num_time_steps 4]))
        regions_infant_coverage_1980_2100_mat(ii,:,:) = sum(countries_infant_coverage_1980_2100_mat .* region_index_mat,1); % for each scenario, sum over all countries in region
        regions_BD_coverage_1980_2100_mat(ii,:,:) = sum(countries_BD_coverage_1980_2100_mat .* region_index_mat,1); % for each scenario, sum over all countries in region
        regions_Tot_Pop_2025_mat(ii,:,:) = sum(countries_Tot_Pop_2025_mat .* region_index_mat,1);
    end
    assert(isequal(size(regions_infant_coverage_1980_2100_mat),[num_regions_global num_time_steps 4]))
    assert(isequal(size(regions_BD_coverage_1980_2100_mat),[num_regions_global num_time_steps 4]))
    assert(isequal(size(regions_Tot_Pop_2025_mat),[num_regions_global num_time_steps 4]))
    assert(all(all(all(regions_infant_coverage_1980_2100_mat>=0))))
    assert(all(all(all(regions_BD_coverage_1980_2100_mat>=0))))
    assert(all(all(all(regions_Tot_Pop_2025_mat>=0))))
    regions_infant_coverage_1980_2100_mat = regions_infant_coverage_1980_2100_mat ./ regions_Tot_Pop_2025_mat;
    regions_BD_coverage_1980_2100_mat = regions_BD_coverage_1980_2100_mat ./ regions_Tot_Pop_2025_mat;
    assert(isequal(size(regions_infant_coverage_1980_2100_mat),[num_regions_global num_time_steps 4]))
    assert(all(all(all(regions_infant_coverage_1980_2100_mat>=0))))
    assert(all(all(all(regions_infant_coverage_1980_2100_mat<=100))))
    assert(isequal(size(regions_BD_coverage_1980_2100_mat),[num_regions_global num_time_steps 4]))
    assert(all(all(all(regions_BD_coverage_1980_2100_mat>=0))))
    assert(all(all(all(regions_BD_coverage_1980_2100_mat<=100))))


    f = figure('Units','centimeters');
    f.Position(3:4) = [17.5 8];
    tiledlayout(1, 2)

    x_mat = repmat(years_vec_01yr',1,4);
    assert(isequal(size(x_mat),[num_time_steps 4]))
    scenarios_of_interest = [...
        status_quo_infant_status_quo_BD_num ...
        status_quo_infant_BD_to_90_num ...
        status_quo_infant_BD_drop_20_num ...
        BD_delayed_slowed_2025_2040_num...
        ];
    colour_mat = scenario_colour_mat(scenarios_of_interest,:);
    num_colours = size(colour_mat,1);
    linestyle_cell_array = scenario_line_style_array(scenarios_of_interest);
    num_linestyles = length(linestyle_cell_array);
    series_index_vec = 1:(num_colours*num_linestyles);


    subplot_1 = nexttile;
    set(gca, 'ColorOrder', colour_mat, 'LineStyleOrder', linestyle_cell_array)

    hold on
    y_mat = squeeze(regions_infant_coverage_1980_2100_mat(end,:,:));
    assert(isequal(size(y_mat),[num_time_steps 4]))
    for jj=1:num_colours
        plot(x_mat(:,jj),y_mat(:,jj),'LineWidth',1.5,'SeriesIndex',series_index_vec((num_colours+1)*jj - num_colours))
    end
    xlim([1980 2100])
    ylim([0 100])
    yticks(0:10:100)
    xlabel('Year')
    ylabel('Coverage (%)')
    title('HepB3')
    annotate_subplot(subplot_1.Position,'e',2)


    subplot_2 = nexttile;
    set(gca, 'ColorOrder', colour_mat, 'LineStyleOrder', linestyle_cell_array)

    hold on
    y_mat = squeeze(regions_BD_coverage_1980_2100_mat(end,:,:));
    assert(isequal(size(y_mat),[num_time_steps 4]))
    for jj=1:num_colours
        plot(x_mat(:,jj),y_mat(:,jj),'LineWidth',1.5,'SeriesIndex',series_index_vec((num_colours+1)*jj - num_colours))
    end
    xlim([1980 2100])
    ylim([0 100])
    yticks(0:10:100)
    xlabel('Year')
    ylabel('Coverage (%)')
    title('Timely HepB-BD vaccination')
    annotate_subplot(subplot_2.Position,'f',2)


    lg = legend(nexttile(1),ListOfScenarioLabels(scenarios_of_interest),'Location','southoutside');
    lg.NumColumns = 2;
    lg.Layout.Tile = 'south';

    assert(strcmp(get(gca, 'FontName'),'Helvetica'))
    subplot_vec = [subplot_1 subplot_2];
    set(subplot_vec,'FontSize',7)



    f = figure('Units','centimeters');
    f.Position(3:4) = [17.5 8];
    tiledlayout(1, 2)

    x_mat = repmat(calandar_years_vec',1,num_regions);
    assert(isequal(size(x_mat),[num_calandar_years num_regions]))
    scenarios_of_interest = [...
        status_quo_infant_status_quo_BD_num ...
        ];


    subplot_1 = nexttile;
    set(gca, 'ColorOrder', region_colour_mat)

    hold on
    y_mat = squeeze(regions_infant_coverage_1980_2100_mat(1:(end-1),1:num_year_divisions:end,1))';
    assert(isequal(size(y_mat),[num_calandar_years num_regions]))
    plot(x_mat,y_mat,'LineWidth',1.5)
    xlim([1980 2100])
    ylim([0 100])
    yticks(0:10:100)
    xlabel('Year')
    ylabel('Coverage (%)')
    title('HepB3')
    annotate_subplot(subplot_1.Position,'a',2)


    subplot_2 = nexttile;
    set(gca, 'ColorOrder', region_colour_mat)

    hold on
    y_mat = squeeze(regions_BD_coverage_1980_2100_mat(1:(end-1),1:num_year_divisions:end,1))';
    assert(isequal(size(y_mat),[num_calandar_years num_regions]))
    plot(x_mat,y_mat,'LineWidth',1.5)
    xlim([1980 2100])
    ylim([0 100])
    yticks(0:10:100)
    xlabel('Year')
    ylabel('Coverage (%)')
    title('Timely HepB-BD vaccination')
    annotate_subplot(subplot_2.Position,'b',2)


    lg = legend(nexttile(1),regions_list,'Location','southoutside','Orientation','horizontal');
    lg.Layout.Tile = 'south';

    assert(strcmp(get(gca, 'FontName'),'Helvetica'))
    subplot_vec = [subplot_1 subplot_2];
    set(subplot_vec,'FontSize',7)


end % end of draw_Figs_1e_1f_Suppl_Figs_2a_2b_fun function





function draw_Figs_2a_2b_2c_2d_fun(...
    num_stochas_runs,...
    calandar_years_vec,...
    ListOfScenarios,ListOfScenarioLabels,scenario_nums_struct,scenario_line_style_array,scenario_colour_mat,...
    countries_list,...
    regions_cell_array)


    years_vec = 2010:2100;
    num_years = length(years_vec);
    ilower = find(calandar_years_vec>=years_vec(1), 1);
    iupper = find(calandar_years_vec>=years_vec(end), 1);
    years_pos_vec = ilower:iupper;


    scenarios_of_interest = [...
        scenario_nums_struct.status_quo_infant_status_quo_BD_num ...
        scenario_nums_struct.status_quo_infant_BD_to_90_num...
        ];


    f = figure('Units','centimeters');
    f.Position(3:4) = [17.5 14];
    tiledlayout(2, 2)

    alphabet_char = ('a':'z')';                 % 26 letter long character
    charlbl = num2cell(alphabet_char(1:4))';    % {'a','b','c','d'}


    subplot_1 = nexttile;
    hold on
    the_denominator = 1e6;
    for scenario_num=scenarios_of_interest

        scenario = ListOfScenarios{scenario_num};


        results_mat = get_model_output(countries_list,years_vec,years_pos_vec,...
            'regions','Global',scenario,scenario_num,regions_cell_array,num_stochas_runs,'Incid_chronic_all_1yr_approx');
        results_mat = results_mat / the_denominator;
        assert(isequal(size(results_mat),[num_stochas_runs num_years]))
        mean_vec = mean(results_mat,1);
        lower_vec = prctile(results_mat,2.5,1);
        upper_vec = prctile(results_mat,97.5,1);
        plotHandles(scenario_num) = plot(years_vec,mean_vec,...
            'Color',scenario_colour_mat(scenario_num,:),'LineStyle',scenario_line_style_array{scenario_num},'LineWidth',1.5);
        plotLabels{scenario_num} = ListOfScenarioLabels{scenario_num};
        plot(years_vec,lower_vec,...
            'Color',scenario_colour_mat(scenario_num,:),'LineStyle','--','LineWidth',1)
        plot(years_vec,upper_vec,...
            'Color',scenario_colour_mat(scenario_num,:),'LineStyle','--','LineWidth',1)
        fill_handle = fill([years_vec fliplr(years_vec)],[upper_vec fliplr(lower_vec)],scenario_colour_mat(scenario_num,:));
        alpha(fill_handle,0.25)

    end
    xlim([years_vec(1) years_vec(end)])
    ylim([0 3e6]/the_denominator)
    xlabel('Year')
    ylabel('Incident cases (millions)')
    annotate_subplot(subplot_1.Position,charlbl{1},2)


    subplot_2 = nexttile;
    hold on
    for scenario_num=scenarios_of_interest

        scenario = ListOfScenarios{scenario_num};


        results_mat = get_model_output(countries_list,years_vec,years_pos_vec,...
            'regions','Global',scenario,scenario_num,regions_cell_array,num_stochas_runs,'NumSAg_chronic_1yr_5_year_olds','Tot_Pop_1yr_5_year_olds');
        results_mat = results_mat * 100;
        assert(isequal(size(results_mat),[num_stochas_runs num_years]))
        mean_vec = mean(results_mat,1);
        lower_vec = prctile(results_mat,2.5,1);
        upper_vec = prctile(results_mat,97.5,1);
        plotHandles(scenario_num) = plot(years_vec,mean_vec,...
            'Color',scenario_colour_mat(scenario_num,:),'LineStyle',scenario_line_style_array{scenario_num},'LineWidth',1.5);
        plotLabels{scenario_num} = ListOfScenarioLabels{scenario_num};
        plot(years_vec,lower_vec,...
            'Color',scenario_colour_mat(scenario_num,:),'LineStyle','--','LineWidth',1)
        plot(years_vec,upper_vec,...
            'Color',scenario_colour_mat(scenario_num,:),'LineStyle','--','LineWidth',1)
        fill_handle = fill([years_vec fliplr(years_vec)],[upper_vec fliplr(lower_vec)],scenario_colour_mat(scenario_num,:));
        alpha(fill_handle,0.25)

    end
    xlim([years_vec(1) years_vec(end)])
    ylim([0 3])
    line([years_vec(1) years_vec(end)],[0.1 0.1],'LineStyle',':','LineWidth',1.5,'Color','k')
    xlabel('Year')
    ylabel('Prevalence (%)')
    annotate_subplot(subplot_2.Position,charlbl{2},2)


    subplot_3 = nexttile;
    hold on
    the_denominator = 1e6;
    for scenario_num=scenarios_of_interest

        scenario = ListOfScenarios{scenario_num};


        results_mat = get_model_output(countries_list,years_vec,years_pos_vec,...
            'regions','Global',scenario,scenario_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx');
        results_mat = results_mat / the_denominator;
        assert(isequal(size(results_mat),[num_stochas_runs num_years]))
        mean_vec = mean(results_mat,1);
        lower_vec = prctile(results_mat,2.5,1);
        upper_vec = prctile(results_mat,97.5,1);
        plotHandles(scenario_num) = plot(years_vec,mean_vec,...
            'Color',scenario_colour_mat(scenario_num,:),'LineStyle',scenario_line_style_array{scenario_num},'LineWidth',1.5);
        plotLabels{scenario_num} = ListOfScenarioLabels{scenario_num};
        plot(years_vec,lower_vec,...
            'Color',scenario_colour_mat(scenario_num,:),'LineStyle','--','LineWidth',1)
        plot(years_vec,upper_vec,...
            'Color',scenario_colour_mat(scenario_num,:),'LineStyle','--','LineWidth',1)
        fill_handle = fill([years_vec fliplr(years_vec)],[upper_vec fliplr(lower_vec)],scenario_colour_mat(scenario_num,:));
        alpha(fill_handle,0.25)

    end
    xlim([years_vec(1) years_vec(end)])
    ylim([0 1e6]/the_denominator)
    xlabel('Year')
    ylabel('HBV-related deaths (millions)')
    annotate_subplot(subplot_3.Position,charlbl{3},2)


    subplot_4 = nexttile;
    hold on
    the_denominator = 1e6;
    for scenario_num=scenarios_of_interest

        scenario = ListOfScenarios{scenario_num};


        results_mat = get_model_output(countries_list,years_vec,years_pos_vec,...
            'regions','Global',scenario,scenario_num,regions_cell_array,num_stochas_runs,'DALYPerYear');
        results_mat = results_mat / the_denominator;
        assert(isequal(size(results_mat),[num_stochas_runs num_years]))
        mean_vec = mean(results_mat,1);
        lower_vec = prctile(results_mat,2.5,1);
        upper_vec = prctile(results_mat,97.5,1);
        plotHandles(scenario_num) = plot(years_vec,mean_vec,...
            'Color',scenario_colour_mat(scenario_num,:),'LineStyle',scenario_line_style_array{scenario_num},'LineWidth',1.5);
        plotLabels{scenario_num} = ListOfScenarioLabels{scenario_num};
        plot(years_vec,lower_vec,...
            'Color',scenario_colour_mat(scenario_num,:),'LineStyle','--','LineWidth',1)
        plot(years_vec,upper_vec,...
            'Color',scenario_colour_mat(scenario_num,:),'LineStyle','--','LineWidth',1)
        fill_handle = fill([years_vec fliplr(years_vec)],[upper_vec fliplr(lower_vec)],scenario_colour_mat(scenario_num,:));
        alpha(fill_handle,0.25)

    end
    xlim([years_vec(1) years_vec(end)])
    ylim([0 25e6]/the_denominator)
    xlabel('Year')
    ylabel('DALYs (millions)')
    annotate_subplot(subplot_4.Position,charlbl{4},2)


    lg = legend(nexttile(3),plotHandles(scenarios_of_interest), plotLabels(scenarios_of_interest),'Orientation','horizontal','Location','southoutside');
    lg.Layout.Tile = 'south';


    assert(strcmp(get(gca, 'FontName'),'Helvetica'))
    subplot_vec = [subplot_1 subplot_2 subplot_3 subplot_4];
    set(subplot_vec,'FontSize',7)


end % end of draw_Figs_2a_2b_2c_2d_fun function






function draw_Fig_3a_fun(...
    num_stochas_runs,...
    cohort_years_vec,num_cohort_years,cohort_years_pos_vec,...
    scenario_nums_struct,...
    num_regions,regions_list,region_colour_mat,...
    countries_list,...
    regions_cell_array)


    status_quo_infant_status_quo_BD_num = scenario_nums_struct.status_quo_infant_status_quo_BD_num;
    status_quo_infant_BD_to_90_num = scenario_nums_struct.status_quo_infant_BD_to_90_num;


    i2020 = find(cohort_years_vec>=2020, 1);
    i2030 = find(cohort_years_vec>=2030, 1);

    [region_status_quo_deaths_mat_mean, ...
        region_status_quo_deaths_mat_lower_2_5, ...
        region_status_quo_deaths_mat_upper_97_5, ... 
        region_BD_90_deaths_mat_mean, ...
        region_BD_90_deaths_mat_lower_2_5, ...
        region_BD_90_deaths_mat_upper_97_5, ... 
        region_excess_deaths_mat_mean, ...
        region_excess_deaths_mat_lower_2_5, ...
        region_excess_deaths_mat_upper_97_5] ... 
        = deal(-99 * ones(num_regions,num_cohort_years));

    region_excess_deaths_mat_lower_mean_upper = -99 * ones(1,3);


    for region_num=1:num_regions
        
        region = regions_list{region_num};

        [region_excess_deaths_mat_mean,region_excess_deaths_mat_lower_2_5,region_excess_deaths_mat_upper_97_5,...
            region_status_quo_deaths_mat_mean,region_status_quo_deaths_mat_lower_2_5,region_status_quo_deaths_mat_upper_97_5,...
            region_BD_90_deaths_mat_mean,region_BD_90_deaths_mat_lower_2_5,region_BD_90_deaths_mat_upper_97_5] = ...
            difference_in_deaths_fun(false,...
            countries_list,cohort_years_vec,num_cohort_years,cohort_years_pos_vec,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD','Status quo infant & BD expansion to 90%',...
            status_quo_infant_status_quo_BD_num,status_quo_infant_BD_to_90_num,...
            region_excess_deaths_mat_mean,region_excess_deaths_mat_lower_2_5,region_excess_deaths_mat_upper_97_5,...
            region_status_quo_deaths_mat_mean,region_status_quo_deaths_mat_lower_2_5,region_status_quo_deaths_mat_upper_97_5,...
            region_BD_90_deaths_mat_mean,region_BD_90_deaths_mat_lower_2_5,region_BD_90_deaths_mat_upper_97_5);

    end % end of region_num for loop


    status_quo_results_mat = get_model_output(countries_list,cohort_years_vec,cohort_years_pos_vec,...
        'regions','Global','Status quo infant & BD',status_quo_infant_status_quo_BD_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts');
    status_quo_results_mat = status_quo_results_mat(:,i2020:i2030);
    assert(isequal(size(status_quo_results_mat),[num_stochas_runs 2030-2020+1]))

    BD_90_results_mat = get_model_output(countries_list,cohort_years_vec,cohort_years_pos_vec,...
        'regions','Global','Status quo infant & BD expansion to 90%',status_quo_infant_BD_to_90_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts');
    BD_90_results_mat = BD_90_results_mat(:,i2020:i2030);
    assert(isequal(size(BD_90_results_mat),[num_stochas_runs 2030-2020+1]))

    diff_mat = status_quo_results_mat - BD_90_results_mat;
    assert(isequal(size(diff_mat),[num_stochas_runs 2030-2020+1]))
    diff_vec = sum(diff_mat,2);
    assert(isequal(size(diff_vec),[num_stochas_runs 1]))
    region_excess_deaths_mat_lower_mean_upper(2) = mean(diff_vec,1);
    region_excess_deaths_mat_lower_mean_upper(1) = prctile(diff_vec,2.5,1);
    region_excess_deaths_mat_lower_mean_upper(3) = prctile(diff_vec,97.5,1);


    assert(isequal(size(region_status_quo_deaths_mat_mean),[num_regions num_cohort_years]))
    assert(isequal(size(region_status_quo_deaths_mat_lower_2_5),[num_regions num_cohort_years]))
    assert(isequal(size(region_status_quo_deaths_mat_upper_97_5),[num_regions num_cohort_years]))
    assert(all(all(region_status_quo_deaths_mat_mean>=0)))
    assert(all(all(region_status_quo_deaths_mat_lower_2_5>=0)))
    assert(all(all(region_status_quo_deaths_mat_upper_97_5>=0)))
    assert(isequal(size(region_BD_90_deaths_mat_mean),[num_regions num_cohort_years]))
    assert(isequal(size(region_BD_90_deaths_mat_lower_2_5),[num_regions num_cohort_years]))
    assert(isequal(size(region_BD_90_deaths_mat_upper_97_5),[num_regions num_cohort_years]))
    assert(all(all(region_BD_90_deaths_mat_mean>=0)))
    assert(all(all(region_BD_90_deaths_mat_lower_2_5>=0)))
    assert(all(all(region_BD_90_deaths_mat_upper_97_5>=0)))
    assert(isequal(size(region_excess_deaths_mat_mean),[num_regions num_cohort_years]))
    assert(isequal(size(region_excess_deaths_mat_lower_2_5),[num_regions num_cohort_years]))
    assert(isequal(size(region_excess_deaths_mat_upper_97_5),[num_regions num_cohort_years]))
    assert(all(all(region_excess_deaths_mat_mean>=0)))
    assert(all(all(region_excess_deaths_mat_lower_2_5>=0)))
    assert(all(all(region_excess_deaths_mat_upper_97_5>=0)))
    assert(isequal(size(region_excess_deaths_mat_lower_mean_upper),[1 3]))
    assert(all(region_excess_deaths_mat_lower_mean_upper>=0))


    figure

    the_denominator = 1e3;
    ylabel_str = 'thousands';

    ax = axes;
    ax.ColorOrder = region_colour_mat;
    hold on
    barplot = bar(cohort_years_vec,region_excess_deaths_mat_mean'/the_denominator,'stacked');
    error_mid_mat = cumsum(region_excess_deaths_mat_mean',2)/the_denominator;
    error_lower_mat = (region_excess_deaths_mat_mean - region_excess_deaths_mat_lower_2_5)'/the_denominator;
    error_upper_mat = (region_excess_deaths_mat_upper_97_5 - region_excess_deaths_mat_mean)'/the_denominator;
    for region_num=1:num_regions
        errorbar(cohort_years_vec,error_mid_mat(:,region_num),error_lower_mat(:,region_num),error_upper_mat(:,region_num),'.k')
    end
    xlim_vec = xlim;
    xlim([cohort_years_vec(1) xlim_vec(2)])
    ylim([0 140e3]/the_denominator)
    xlabel('Birth cohort')
    ylabel(['HBV-related deaths averted (' ylabel_str ')'])
    lgd = legend(regions_list,'Location','northeast');
    %title({'Additional HBV Deaths in status quo';'Compared to timely HepB-BD expansion to 90%'})
    annotate_subplot(get(gca,'Position'),'a',1)


    assert(strcmp(get(gca, 'FontName'),'Helvetica'))
    set(gca,'FontSize',7)


    disp('Deaths averted globally in scenario 2 relative to scenario 1 in the 2020 and 2030 birth cohorts')
    mean_num = region_excess_deaths_mat_lower_mean_upper(2);
    lower_num = region_excess_deaths_mat_lower_mean_upper(1);
    upper_num = region_excess_deaths_mat_lower_mean_upper(3);
    assert(mean_num>lower_num)
    assert(mean_num<upper_num)
    disp([num2str(round(mean_num)) ' (' num2str(round(lower_num)) ', ' num2str(round(upper_num)) ')'])


end % end of draw_Fig_3a_fun function






function [countries_list, perc_global_total_excess_deaths_vec, perc_global_reduc_in_deaths_vec] = draw_Fig_3b_fun(...
    num_stochas_runs,...
    cohort_years_vec,cohort_years_pos_vec,...
    scenario_nums_struct,...
    num_countries,countries_list,...
    countries_cell_array)


    status_quo_infant_status_quo_BD_num = scenario_nums_struct.status_quo_infant_status_quo_BD_num;
    status_quo_infant_BD_to_90_num = scenario_nums_struct.status_quo_infant_BD_to_90_num;


    cohort_years_vec_1yr = 2020:2030;
    num_years = length(cohort_years_vec_1yr);
    ilower = find(cohort_years_vec>=cohort_years_vec_1yr(1), 1);
    iupper = find(cohort_years_vec>=cohort_years_vec_1yr(end), 1);
    cohort_years_pos_vec_1yr = (cohort_years_pos_vec(ilower)):(cohort_years_pos_vec(iupper));

    [country_status_quo_deaths_mat_mean, ...
        country_status_quo_deaths_mat_lower_2_5, ...
        country_status_quo_deaths_mat_upper_97_5, ...
        country_BD_90_deaths_mat_mean, ...
        country_BD_90_deaths_mat_lower_2_5, ...
        country_BD_90_deaths_mat_upper_97_5, ...
        country_excess_deaths_mat_mean, ...
        country_excess_deaths_mat_lower_2_5, ...
        country_excess_deaths_mat_upper_97_5] = ...
        deal(-99 * ones(num_countries,num_years));


    for country_num=1:num_countries
        
        ISO = countries_list{country_num};

        [country_excess_deaths_mat_mean,country_excess_deaths_mat_lower_2_5,country_excess_deaths_mat_upper_97_5,...
            country_status_quo_deaths_mat_mean,country_status_quo_deaths_mat_lower_2_5,country_status_quo_deaths_mat_upper_97_5,...
            country_BD_90_deaths_mat_mean,country_BD_90_deaths_mat_lower_2_5,country_BD_90_deaths_mat_upper_97_5] = ...
            difference_in_deaths_fun(false,...
            countries_list,cohort_years_vec_1yr,num_years,cohort_years_pos_vec_1yr,'countries',ISO,country_num,countries_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD','Status quo infant & BD expansion to 90%',...
            status_quo_infant_status_quo_BD_num,status_quo_infant_BD_to_90_num,...
            country_excess_deaths_mat_mean,country_excess_deaths_mat_lower_2_5,country_excess_deaths_mat_upper_97_5,...
            country_status_quo_deaths_mat_mean,country_status_quo_deaths_mat_lower_2_5,country_status_quo_deaths_mat_upper_97_5,...
            country_BD_90_deaths_mat_mean,country_BD_90_deaths_mat_lower_2_5,country_BD_90_deaths_mat_upper_97_5);

    end % end of country_num for loop


    assert(isequal(size(country_status_quo_deaths_mat_mean),[num_countries num_years]))
    assert(isequal(size(country_status_quo_deaths_mat_lower_2_5),[num_countries num_years]))
    assert(isequal(size(country_status_quo_deaths_mat_upper_97_5),[num_countries num_years]))
    assert(all(all(country_status_quo_deaths_mat_mean>=0)))
    assert(all(all(country_status_quo_deaths_mat_lower_2_5>=0)))
    assert(all(all(country_status_quo_deaths_mat_upper_97_5>=0)))
    assert(isequal(size(country_BD_90_deaths_mat_mean),[num_countries num_years]))
    assert(isequal(size(country_BD_90_deaths_mat_lower_2_5),[num_countries num_years]))
    assert(isequal(size(country_BD_90_deaths_mat_upper_97_5),[num_countries num_years]))
    assert(all(all(country_BD_90_deaths_mat_mean>=0)))
    assert(all(all(country_BD_90_deaths_mat_lower_2_5>=0)))
    assert(all(all(country_BD_90_deaths_mat_upper_97_5>=0)))
    assert(isequal(size(country_excess_deaths_mat_mean),[num_countries num_years]))
    assert(isequal(size(country_excess_deaths_mat_lower_2_5),[num_countries num_years]))
    assert(isequal(size(country_excess_deaths_mat_upper_97_5),[num_countries num_years]))
    assert(all(all(country_excess_deaths_mat_mean>=0)))
    assert(all(all(country_excess_deaths_mat_lower_2_5>=0)))
    assert(all(all(country_excess_deaths_mat_upper_97_5>=0)))

    country_status_quo_deaths_mat = country_status_quo_deaths_mat_mean;
    country_BD_90_deaths_mat = country_BD_90_deaths_mat_mean;
    country_excess_deaths_mat = country_excess_deaths_mat_mean;

    sum_mean_num_status_quo_deaths_vec = sum(country_status_quo_deaths_mat,2);
    sum_mean_num_BD_90_deaths_vec = sum(country_BD_90_deaths_mat,2);
    sum_mean_num_excess_deaths_vec = sum(country_excess_deaths_mat,2);
    assert(length(sum_mean_num_status_quo_deaths_vec)==num_countries)
    assert(length(sum_mean_num_BD_90_deaths_vec)==num_countries)
    assert(length(sum_mean_num_excess_deaths_vec)==num_countries)
    %[sorted_vec pos_in_sorted_vec] = sort(sum_mean_num_excess_deaths_vec);
    perc_global_total_excess_deaths_vec = sum_mean_num_excess_deaths_vec ./ sum(sum_mean_num_excess_deaths_vec) * 100;
    disp("Use the cell array countries_list and vector perc_global_total_excess_deaths_vec to draw a map in Excel for Fig. 3b.")
    perc_global_reduc_in_deaths_vec = (sum_mean_num_status_quo_deaths_vec - sum_mean_num_BD_90_deaths_vec) ./ sum_mean_num_status_quo_deaths_vec * 100;
    disp("Use the cell array countries_list and vector perc_global_reduc_in_deaths_vec to draw a map in Excel for Supplementary Fig. 3.")


end % end of draw_Fig_3b_fun function






function fig_struct = prepare_Fig_4a_data_struc_fun(...
    num_stochas_runs,...
    calandar_years_vec,calandar_years_pos_vec,...
    ListOfScenarioLabels,scenario_nums_struct,scenario_colour_mat,...
    num_regions,regions_list,num_regions_global,regions_global_list,index_global,...
    countries_list,...
    regions_cell_array)


    status_quo_infant_status_quo_BD_num = scenario_nums_struct.status_quo_infant_status_quo_BD_num;
    status_quo_infant_BD_to_90_num = scenario_nums_struct.status_quo_infant_BD_to_90_num;


    [regions_year_of_elim_prev_status_quo_mat_lower_median_upper,...
        regions_year_of_elim_prev_BD_90_mat_lower_median_upper] = deal(zeros(num_regions_global,3));


    for region_num=1:num_regions_global
        
        region = regions_global_list{region_num};

        regions_year_of_elim_prev_status_quo_mat_lower_median_upper = ...
            find_year_of_elim_fun(countries_list,calandar_years_vec,calandar_years_pos_vec,...
            'regions',region,region_num,regions_cell_array,num_stochas_runs,0.1,...
            'Status quo infant & BD',status_quo_infant_status_quo_BD_num,...
            regions_year_of_elim_prev_status_quo_mat_lower_median_upper);

        regions_year_of_elim_prev_BD_90_mat_lower_median_upper = ...
            find_year_of_elim_fun(countries_list,calandar_years_vec,calandar_years_pos_vec,...
            'regions',region,region_num,regions_cell_array,num_stochas_runs,0.1,...
            'Status quo infant & BD expansion to 90%',status_quo_infant_BD_to_90_num,...
            regions_year_of_elim_prev_BD_90_mat_lower_median_upper);

    end % end of region_num for loop


    assert(isequal(size(regions_year_of_elim_prev_status_quo_mat_lower_median_upper),[num_regions_global 3]))
    assert(all(all(regions_year_of_elim_prev_status_quo_mat_lower_median_upper>=0)))
    assert(isequal(size(regions_year_of_elim_prev_BD_90_mat_lower_median_upper),[num_regions_global 3]))
    assert(all(all(regions_year_of_elim_prev_BD_90_mat_lower_median_upper>=0)))


    row_indices_not_global = 1:num_regions_global~=index_global;


    fig_struct.ListOfScenarioLabels = ListOfScenarioLabels;
    fig_struct.scenario_colour_mat = scenario_colour_mat;
    fig_struct.status_quo_infant_status_quo_BD_num = status_quo_infant_status_quo_BD_num;
    fig_struct.status_quo_infant_BD_to_90_num = status_quo_infant_BD_to_90_num;
    fig_struct.regions_list = regions_list;
    fig_struct.num_regions = num_regions;
    fig_struct.row_indices_not_global = row_indices_not_global;
    fig_struct.regions_year_of_elim_prev_status_quo_mat_lower_median_upper = regions_year_of_elim_prev_status_quo_mat_lower_median_upper;
    fig_struct.regions_year_of_elim_prev_BD_90_mat_lower_median_upper = regions_year_of_elim_prev_BD_90_mat_lower_median_upper;


    disp('Global year of elimination (HBsAg prevalence falls below 0.1% in five-year-olds) in scenario 2')
    median_num = round(regions_year_of_elim_prev_BD_90_mat_lower_median_upper(index_global,2));
    lower_num = round(regions_year_of_elim_prev_BD_90_mat_lower_median_upper(index_global,1));
    upper_num = round(regions_year_of_elim_prev_BD_90_mat_lower_median_upper(index_global,3));
    adjust_years_fun(lower_num,median_num,upper_num)


end % end of prepare_Fig_4a_data_struc_fun function





function fig_struct = prepare_Fig_4b_data_struc_fun(...
    num_stochas_runs,num_year_divisions,...
    calandar_years_vec,num_calandar_years,calandar_years_pos_vec,...
    ListOfScenarios,scenario_nums_struct,...
    num_regions,regions_list,Regions_map,region_colour_mat,...
    num_countries,countries_list,...
    countries_cell_array,vaccination_map_default)


    status_quo_infant_status_quo_BD_num = scenario_nums_struct.status_quo_infant_status_quo_BD_num;
    status_quo_infant_BD_to_25_num = scenario_nums_struct.status_quo_infant_BD_to_25_num;
    status_quo_infant_BD_to_50_num = scenario_nums_struct.status_quo_infant_BD_to_50_num;
    status_quo_infant_BD_to_75_num = scenario_nums_struct.status_quo_infant_BD_to_75_num;
    status_quo_infant_BD_to_90_num = scenario_nums_struct.status_quo_infant_BD_to_90_num;


    i2020 = find(calandar_years_vec>=2020, 1);
    i2030 = find(calandar_years_vec>=2030, 1);

    [countries_BD_coverage_2020_2030_status_quo_mat, ...
        countries_BD_coverage_2020_2030_BD_25_mat, ...
        countries_BD_coverage_2020_2030_BD_50_mat, ...
        countries_BD_coverage_2020_2030_BD_75_mat, ...
        countries_BD_coverage_2020_2030_BD_90_mat] = ...
        deal(zeros(num_countries,2));
    [countries_Tot_Pop_2025_status_quo_mat, ...
        countries_Tot_Pop_2025_BD_25_mat, ...
        countries_Tot_Pop_2025_BD_50_mat, ...
        countries_Tot_Pop_2025_BD_75_mat, ...
        countries_Tot_Pop_2025_BD_90_mat] = ...
        deal(zeros(num_countries,1));

    [countries_prev_chronic_cases_5_year_olds_status_quo_mat, ...
        countries_prev_chronic_cases_5_year_olds_BD_25_mat, ...
        countries_prev_chronic_cases_5_year_olds_BD_50_mat, ...
        countries_prev_chronic_cases_5_year_olds_BD_75_mat, ...
        countries_prev_chronic_cases_5_year_olds_BD_90_mat] = ...
        deal(-99 * ones(num_countries,num_stochas_runs,num_calandar_years));
    [countries_Tot_Pop_5_year_olds_status_quo_mat, ...
        countries_Tot_Pop_5_year_olds_BD_25_mat, ...
        countries_Tot_Pop_5_year_olds_BD_50_mat, ...
        countries_Tot_Pop_5_year_olds_BD_75_mat, ...
        countries_Tot_Pop_5_year_olds_BD_90_mat] = ...
        deal(-99 * ones(num_countries,num_stochas_runs,num_calandar_years));


    for country_num=1:num_countries

        ISO = countries_list{country_num};

        [countries_BD_coverage_2020_2030_status_quo_mat,countries_Tot_Pop_2025_status_quo_mat,...
            countries_prev_chronic_cases_5_year_olds_status_quo_mat,countries_Tot_Pop_5_year_olds_status_quo_mat] = ...
            get_BD_NumSAg_Tot_Pop(num_countries,countries_list,...
            num_year_divisions,calandar_years_vec,num_calandar_years,calandar_years_pos_vec,i2020,i2030,...
            'countries',ISO,country_num,countries_cell_array,vaccination_map_default,num_stochas_runs,...
            ListOfScenarios,'Status quo infant & BD',status_quo_infant_status_quo_BD_num,...
            countries_BD_coverage_2020_2030_status_quo_mat,countries_Tot_Pop_2025_status_quo_mat,...
            countries_prev_chronic_cases_5_year_olds_status_quo_mat,countries_Tot_Pop_5_year_olds_status_quo_mat);

        [countries_BD_coverage_2020_2030_BD_25_mat,countries_Tot_Pop_2025_BD_25_mat,...
            countries_prev_chronic_cases_5_year_olds_BD_25_mat,countries_Tot_Pop_5_year_olds_BD_25_mat] = ...
            get_BD_NumSAg_Tot_Pop(num_countries,countries_list,...
            num_year_divisions,calandar_years_vec,num_calandar_years,calandar_years_pos_vec,i2020,i2030,...
            'countries',ISO,country_num,countries_cell_array,vaccination_map_default,num_stochas_runs,...
            ListOfScenarios,'Status quo infant & BD expansion to 25%',status_quo_infant_BD_to_25_num,...
            countries_BD_coverage_2020_2030_BD_25_mat,countries_Tot_Pop_2025_BD_25_mat,...
            countries_prev_chronic_cases_5_year_olds_BD_25_mat,countries_Tot_Pop_5_year_olds_BD_25_mat);

        [countries_BD_coverage_2020_2030_BD_50_mat,countries_Tot_Pop_2025_BD_50_mat,...
            countries_prev_chronic_cases_5_year_olds_BD_50_mat,countries_Tot_Pop_5_year_olds_BD_50_mat] = ...
            get_BD_NumSAg_Tot_Pop(num_countries,countries_list,...
            num_year_divisions,calandar_years_vec,num_calandar_years,calandar_years_pos_vec,i2020,i2030,...
            'countries',ISO,country_num,countries_cell_array,vaccination_map_default,num_stochas_runs,...
            ListOfScenarios,'Status quo infant & BD expansion to 50%',status_quo_infant_BD_to_50_num,...
            countries_BD_coverage_2020_2030_BD_50_mat,countries_Tot_Pop_2025_BD_50_mat,...
            countries_prev_chronic_cases_5_year_olds_BD_50_mat,countries_Tot_Pop_5_year_olds_BD_50_mat);

        [countries_BD_coverage_2020_2030_BD_75_mat,countries_Tot_Pop_2025_BD_75_mat,...
            countries_prev_chronic_cases_5_year_olds_BD_75_mat,countries_Tot_Pop_5_year_olds_BD_75_mat] = ...
            get_BD_NumSAg_Tot_Pop(num_countries,countries_list,...
            num_year_divisions,calandar_years_vec,num_calandar_years,calandar_years_pos_vec,i2020,i2030,...
            'countries',ISO,country_num,countries_cell_array,vaccination_map_default,num_stochas_runs,...
            ListOfScenarios,'Status quo infant & BD expansion to 75%',status_quo_infant_BD_to_75_num,...
            countries_BD_coverage_2020_2030_BD_75_mat,countries_Tot_Pop_2025_BD_75_mat,...
            countries_prev_chronic_cases_5_year_olds_BD_75_mat,countries_Tot_Pop_5_year_olds_BD_75_mat);

        [countries_BD_coverage_2020_2030_BD_90_mat,countries_Tot_Pop_2025_BD_90_mat,...
            countries_prev_chronic_cases_5_year_olds_BD_90_mat,countries_Tot_Pop_5_year_olds_BD_90_mat] = ...
            get_BD_NumSAg_Tot_Pop(num_countries,countries_list,...
            num_year_divisions,calandar_years_vec,num_calandar_years,calandar_years_pos_vec,i2020,i2030,...
            'countries',ISO,country_num,countries_cell_array,vaccination_map_default,num_stochas_runs,...
            ListOfScenarios,'Status quo infant & BD expansion to 90%',status_quo_infant_BD_to_90_num,...
            countries_BD_coverage_2020_2030_BD_90_mat,countries_Tot_Pop_2025_BD_90_mat,...
            countries_prev_chronic_cases_5_year_olds_BD_90_mat,countries_Tot_Pop_5_year_olds_BD_90_mat);

    end % end of country_num for loop


    assert(isequal(size(countries_BD_coverage_2020_2030_status_quo_mat),[num_countries 2]))
    assert(isequal(size(countries_BD_coverage_2020_2030_BD_25_mat),[num_countries 2]))
    assert(isequal(size(countries_BD_coverage_2020_2030_BD_50_mat),[num_countries 2]))
    assert(isequal(size(countries_BD_coverage_2020_2030_BD_75_mat),[num_countries 2]))
    assert(isequal(size(countries_BD_coverage_2020_2030_BD_90_mat),[num_countries 2]))
    assert(all(all(countries_BD_coverage_2020_2030_status_quo_mat>=0)))
    assert(all(all(countries_BD_coverage_2020_2030_BD_25_mat>=0)))
    assert(all(all(countries_BD_coverage_2020_2030_BD_50_mat>=0)))
    assert(all(all(countries_BD_coverage_2020_2030_BD_75_mat>=0)))
    assert(all(all(countries_BD_coverage_2020_2030_BD_90_mat>=0)))
    assert(all(all(countries_BD_coverage_2020_2030_status_quo_mat<=100)))
    assert(all(all(countries_BD_coverage_2020_2030_BD_25_mat<=100)))
    assert(all(all(countries_BD_coverage_2020_2030_BD_50_mat<=100)))
    assert(all(all(countries_BD_coverage_2020_2030_BD_75_mat<=100)))
    assert(all(all(countries_BD_coverage_2020_2030_BD_90_mat<=100)))
    assert(all(all(countries_BD_coverage_2020_2030_status_quo_mat<=countries_BD_coverage_2020_2030_BD_25_mat)))
    assert(all(all(countries_BD_coverage_2020_2030_BD_25_mat<=countries_BD_coverage_2020_2030_BD_50_mat)))
    assert(all(all(countries_BD_coverage_2020_2030_BD_50_mat<=countries_BD_coverage_2020_2030_BD_75_mat)))
    assert(all(all(countries_BD_coverage_2020_2030_BD_75_mat<=countries_BD_coverage_2020_2030_BD_90_mat)))
    assert(isequal(size(countries_Tot_Pop_2025_status_quo_mat),[num_countries 1]))
    assert(isequal(size(countries_Tot_Pop_2025_BD_25_mat),[num_countries 1]))
    assert(isequal(size(countries_Tot_Pop_2025_BD_50_mat),[num_countries 1]))
    assert(isequal(size(countries_Tot_Pop_2025_BD_75_mat),[num_countries 1]))
    assert(isequal(size(countries_Tot_Pop_2025_BD_90_mat),[num_countries 1]))
    assert(all(all(countries_Tot_Pop_2025_status_quo_mat>=0)))
    assert(all(all(countries_Tot_Pop_2025_BD_25_mat>=0)))
    assert(all(all(countries_Tot_Pop_2025_BD_50_mat>=0)))
    assert(all(all(countries_Tot_Pop_2025_BD_75_mat>=0)))
    assert(all(all(countries_Tot_Pop_2025_BD_90_mat>=0)))
    assert(isequal(size(countries_prev_chronic_cases_5_year_olds_status_quo_mat),[num_countries num_stochas_runs num_calandar_years]))
    assert(isequal(size(countries_prev_chronic_cases_5_year_olds_BD_25_mat),[num_countries num_stochas_runs num_calandar_years]))
    assert(isequal(size(countries_prev_chronic_cases_5_year_olds_BD_50_mat),[num_countries num_stochas_runs num_calandar_years]))
    assert(isequal(size(countries_prev_chronic_cases_5_year_olds_BD_75_mat),[num_countries num_stochas_runs num_calandar_years]))
    assert(isequal(size(countries_prev_chronic_cases_5_year_olds_BD_90_mat),[num_countries num_stochas_runs num_calandar_years]))
    assert(all(all(all(countries_prev_chronic_cases_5_year_olds_status_quo_mat>=0))))
    assert(all(all(all(countries_prev_chronic_cases_5_year_olds_BD_25_mat>=0))))
    assert(all(all(all(countries_prev_chronic_cases_5_year_olds_BD_50_mat>=0))))
    assert(all(all(all(countries_prev_chronic_cases_5_year_olds_BD_75_mat>=0))))
    assert(all(all(all(countries_prev_chronic_cases_5_year_olds_BD_90_mat>=0))))
    assert(all(all(all(countries_prev_chronic_cases_5_year_olds_status_quo_mat>=countries_prev_chronic_cases_5_year_olds_BD_25_mat))))
    assert(all(all(all(countries_prev_chronic_cases_5_year_olds_BD_25_mat>=countries_prev_chronic_cases_5_year_olds_BD_50_mat))))
    assert(all(all(all(countries_prev_chronic_cases_5_year_olds_BD_50_mat>=countries_prev_chronic_cases_5_year_olds_BD_75_mat))))
    assert(all(all(all(countries_prev_chronic_cases_5_year_olds_BD_75_mat>=countries_prev_chronic_cases_5_year_olds_BD_90_mat))))
    assert(isequal(size(countries_Tot_Pop_5_year_olds_status_quo_mat),[num_countries num_stochas_runs num_calandar_years]))
    assert(isequal(size(countries_Tot_Pop_5_year_olds_BD_25_mat),[num_countries num_stochas_runs num_calandar_years]))
    assert(isequal(size(countries_Tot_Pop_5_year_olds_BD_50_mat),[num_countries num_stochas_runs num_calandar_years]))
    assert(isequal(size(countries_Tot_Pop_5_year_olds_BD_75_mat),[num_countries num_stochas_runs num_calandar_years]))
    assert(isequal(size(countries_Tot_Pop_5_year_olds_BD_90_mat),[num_countries num_stochas_runs num_calandar_years]))
    assert(all(all(all(countries_Tot_Pop_5_year_olds_status_quo_mat>=0))))
    assert(all(all(all(countries_Tot_Pop_5_year_olds_BD_25_mat>=0))))
    assert(all(all(all(countries_Tot_Pop_5_year_olds_BD_50_mat>=0))))
    assert(all(all(all(countries_Tot_Pop_5_year_olds_BD_75_mat>=0))))
    assert(all(all(all(countries_Tot_Pop_5_year_olds_BD_90_mat>=0))))


    countries_prev_chronic_cases_5_year_olds_mat = cat(4,...
        countries_prev_chronic_cases_5_year_olds_BD_25_mat,...
        countries_prev_chronic_cases_5_year_olds_BD_50_mat,...
        countries_prev_chronic_cases_5_year_olds_BD_75_mat,...
        countries_prev_chronic_cases_5_year_olds_BD_90_mat);
    countries_Tot_Pop_5_year_olds_mat = cat(4,...
        countries_Tot_Pop_5_year_olds_BD_25_mat,...
        countries_Tot_Pop_5_year_olds_BD_50_mat,...
        countries_Tot_Pop_5_year_olds_BD_75_mat,...
        countries_Tot_Pop_5_year_olds_BD_90_mat);
    assert(isequal(size(countries_prev_chronic_cases_5_year_olds_mat),[num_countries num_stochas_runs num_calandar_years 4]))
    assert(all(all(all(all(countries_prev_chronic_cases_5_year_olds_mat>=0)))))
    assert(isequal(size(countries_Tot_Pop_5_year_olds_mat),[num_countries num_stochas_runs num_calandar_years 4]))
    assert(all(all(all(all(countries_Tot_Pop_5_year_olds_mat>=0)))))
    countries_BD_coverage_2020_2030_mat = [...
        diff(countries_BD_coverage_2020_2030_BD_25_mat,1,2),...
        diff(countries_BD_coverage_2020_2030_BD_50_mat,1,2),...
        diff(countries_BD_coverage_2020_2030_BD_75_mat,1,2),...
        diff(countries_BD_coverage_2020_2030_BD_90_mat,1,2)]; 
    countries_Tot_Pop_2025_mat = [...
        countries_Tot_Pop_2025_BD_25_mat,...
        countries_Tot_Pop_2025_BD_50_mat,...
        countries_Tot_Pop_2025_BD_75_mat,...
        countries_Tot_Pop_2025_BD_90_mat];
    assert(isequal(size(countries_BD_coverage_2020_2030_mat),[num_countries 4]))
    assert(all(all(countries_BD_coverage_2020_2030_mat>=0)))
    assert(all(all(countries_BD_coverage_2020_2030_mat<=100)))
    assert(isequal(size(countries_Tot_Pop_2025_mat),[num_countries 4]))
    assert(all(all(countries_Tot_Pop_2025_mat>=0)))
    countries_BD_coverage_2020_2030_mat = countries_BD_coverage_2020_2030_mat .* countries_Tot_Pop_2025_mat; % weight timely HepB-BD coverage by population size in 2025
    assert(isequal(size(countries_BD_coverage_2020_2030_mat),[num_countries 4]))

    assert(isequal(size(countries_list),[1 num_countries]))
    countries_list_regions = cellfun(@(x) Regions_map(x),countries_list,'UniformOutput',false);
    assert(isequal(size(countries_list_regions),[1 num_countries]))
    region_index_map = containers.Map(regions_list,{1,2,3,4,5,6});
    countries_index_vec = cellfun(@(x) region_index_map(x),countries_list_regions);
    assert(isequal(size(countries_index_vec),[1 num_countries]))

    regions_year_of_elim_prev_mat = -99 * ones(num_regions,num_stochas_runs,4);
    regions_BD_coverage_2020_2030_mat = -99 * ones(num_regions,4);
    regions_Tot_Pop_2025_mat = -99 * ones(num_regions,4);
    for ii=1:num_regions
        region_index_vec = countries_index_vec == ii;
        assert(isequal(size(region_index_vec),[1 num_countries]))
        assert(all(strcmp(countries_list_regions(region_index_vec),regions_list(ii)))) 

        region_index_mat = repmat(region_index_vec',1,num_stochas_runs,num_calandar_years,4);
        region_prev_chronic_cases_5_year_olds_mat = squeeze(sum(countries_prev_chronic_cases_5_year_olds_mat .* region_index_mat,1));
        region_Tot_Pop_5_year_olds_mat = squeeze(sum(countries_Tot_Pop_5_year_olds_mat .* region_index_mat,1));
        assert(isequal(size(region_prev_chronic_cases_5_year_olds_mat),[num_stochas_runs num_calandar_years 4]))
        assert(isequal(size(region_Tot_Pop_5_year_olds_mat),[num_stochas_runs num_calandar_years 4]))
        region_prev_chronic_5_year_olds_mat = permute(region_prev_chronic_cases_5_year_olds_mat ./ region_Tot_Pop_5_year_olds_mat,[3 1 2]) * 100;
        assert(isequal(size(region_prev_chronic_5_year_olds_mat),[4 num_stochas_runs num_calandar_years]))
        region_year_of_elim_prev_mat = -99 * ones(num_stochas_runs,4);
        for jj=1:num_stochas_runs
            region_year_of_elim_prev_mat(jj,:) = arrayfun(@(xx) ...
                assign_year_of_elim(calandar_years_vec,region_prev_chronic_5_year_olds_mat(xx,jj,:),0.1), ...
                1:size(region_prev_chronic_5_year_olds_mat,1));
        end
        assert(isequal(size(region_year_of_elim_prev_mat),[num_stochas_runs 4]))
        regions_year_of_elim_prev_mat(ii,:,:) = region_year_of_elim_prev_mat;

        region_index_mat = repmat(region_index_vec',1,4);
        regions_BD_coverage_2020_2030_mat(ii,:) = sum(countries_BD_coverage_2020_2030_mat .* region_index_mat,1);
        regions_Tot_Pop_2025_mat(ii,:) = sum(countries_Tot_Pop_2025_mat .* region_index_mat,1);
    end
    assert(isequal(size(regions_BD_coverage_2020_2030_mat),[num_regions 4]))
    assert(isequal(size(regions_Tot_Pop_2025_mat),[num_regions 4]))
    assert(all(all(regions_BD_coverage_2020_2030_mat>=0)))
    assert(all(all(regions_Tot_Pop_2025_mat>=0)))
    regions_BD_coverage_2020_2030_mat = regions_BD_coverage_2020_2030_mat ./ regions_Tot_Pop_2025_mat;
    assert(all(all(regions_BD_coverage_2020_2030_mat>=0)))
    assert(all(all(regions_BD_coverage_2020_2030_mat<=100)))

    assert(isequal(size(regions_year_of_elim_prev_mat),[num_regions num_stochas_runs 4]))
    assert(all(all(all(regions_year_of_elim_prev_mat>=0))))
    median_regions_year_of_elim_prev_mat = round([...
        median(regions_year_of_elim_prev_mat(:,:,1),2),...
        median(regions_year_of_elim_prev_mat(:,:,2),2),...
        median(regions_year_of_elim_prev_mat(:,:,3),2),...
        median(regions_year_of_elim_prev_mat(:,:,4),2)]);
    lower_regions_year_of_elim_prev_mat = round([...
        prctile(regions_year_of_elim_prev_mat(:,:,1),2.5,2),...
        prctile(regions_year_of_elim_prev_mat(:,:,2),2.5,2),...
        prctile(regions_year_of_elim_prev_mat(:,:,3),2.5,2),...
        prctile(regions_year_of_elim_prev_mat(:,:,4),2.5,2)]);
    upper_regions_year_of_elim_prev_mat = round([...
        prctile(regions_year_of_elim_prev_mat(:,:,1),97.5,2),...
        prctile(regions_year_of_elim_prev_mat(:,:,2),97.5,2),...
        prctile(regions_year_of_elim_prev_mat(:,:,3),97.5,2),...
        prctile(regions_year_of_elim_prev_mat(:,:,4),97.5,2)]);
    assert(isequal(size(median_regions_year_of_elim_prev_mat),[num_regions 4]))
    assert(isequal(size(lower_regions_year_of_elim_prev_mat),[num_regions 4]))
    assert(isequal(size(upper_regions_year_of_elim_prev_mat),[num_regions 4]))
    assert(all(all(lower_regions_year_of_elim_prev_mat <= median_regions_year_of_elim_prev_mat)))
    assert(all(all(median_regions_year_of_elim_prev_mat <= upper_regions_year_of_elim_prev_mat)))


    fig_struct.regions_list = regions_list;
    fig_struct.num_regions = num_regions;
    fig_struct.region_colour_mat = region_colour_mat;
    fig_struct.regions_BD_coverage_2020_2030_mat = regions_BD_coverage_2020_2030_mat;
    fig_struct.median_regions_year_of_elim_prev_mat = median_regions_year_of_elim_prev_mat;
    fig_struct.lower_regions_year_of_elim_prev_mat = lower_regions_year_of_elim_prev_mat;
    fig_struct.upper_regions_year_of_elim_prev_mat = upper_regions_year_of_elim_prev_mat;


    index_AFRO = find(strcmp('AFRO',regions_list));
    disp('Year of elimination (HBsAg prevalence falls below 0.1% in five-year-olds) in AFRO in scenario 2a')
    median_num = round(median_regions_year_of_elim_prev_mat(index_AFRO,1));
    lower_num = round(lower_regions_year_of_elim_prev_mat(index_AFRO,1));
    upper_num = round(upper_regions_year_of_elim_prev_mat(index_AFRO,1));
    adjust_years_fun(lower_num,median_num,upper_num)

    disp('Year of elimination (HBsAg prevalence falls below 0.1% in five-year-olds) in AFRO in scenario 2b')
    median_num = round(median_regions_year_of_elim_prev_mat(index_AFRO,2));
    lower_num = round(lower_regions_year_of_elim_prev_mat(index_AFRO,2));
    upper_num = round(upper_regions_year_of_elim_prev_mat(index_AFRO,2));
    adjust_years_fun(lower_num,median_num,upper_num)

    disp('Year of elimination (HBsAg prevalence falls below 0.1% in five-year-olds) in AFRO in scenario 2c')
    median_num = round(median_regions_year_of_elim_prev_mat(index_AFRO,3));
    lower_num = round(lower_regions_year_of_elim_prev_mat(index_AFRO,3));
    upper_num = round(upper_regions_year_of_elim_prev_mat(index_AFRO,3));
    adjust_years_fun(lower_num,median_num,upper_num)

    disp('Year of elimination (HBsAg prevalence falls below 0.1% in five-year-olds) in AFRO in scenario 2d')
    median_num = round(median_regions_year_of_elim_prev_mat(index_AFRO,4));
    lower_num = round(lower_regions_year_of_elim_prev_mat(index_AFRO,4));
    upper_num = round(upper_regions_year_of_elim_prev_mat(index_AFRO,4));
    adjust_years_fun(lower_num,median_num,upper_num)


end % end of prepare_Fig_4b_data_struc_fun function




function draw_Figs_4a_4b_fun(Fig_4a_data_struc,Fig_4b_data_struc)


    f = figure('Units','centimeters');
    f.Position(3:4) = [17.5 8];
    tiledlayout(1, 2)


    subplot_1 = nexttile;
    ListOfScenarioLabels = Fig_4a_data_struc.ListOfScenarioLabels;
    scenario_colour_mat = Fig_4a_data_struc.scenario_colour_mat;
    status_quo_infant_status_quo_BD_num = Fig_4a_data_struc.status_quo_infant_status_quo_BD_num;
    status_quo_infant_BD_to_90_num = Fig_4a_data_struc.status_quo_infant_BD_to_90_num;
    regions_list = Fig_4a_data_struc.regions_list;
    num_regions = Fig_4a_data_struc.num_regions;
    row_indices_not_global = Fig_4a_data_struc.row_indices_not_global;
    regions_year_of_elim_prev_status_quo_mat_lower_median_upper = Fig_4a_data_struc.regions_year_of_elim_prev_status_quo_mat_lower_median_upper;
    regions_year_of_elim_prev_BD_90_mat_lower_median_upper = Fig_4a_data_struc.regions_year_of_elim_prev_BD_90_mat_lower_median_upper;


    hold on
    barplot_x_mat = (1:num_regions)';
    assert(isequal(size(barplot_x_mat),[num_regions 1]))
    barplot = barh(barplot_x_mat,regions_year_of_elim_prev_BD_90_mat_lower_median_upper(row_indices_not_global,2));
    barplot.FaceColor = scenario_colour_mat(status_quo_infant_BD_to_90_num,:);
    error_mid_mat = regions_year_of_elim_prev_BD_90_mat_lower_median_upper(row_indices_not_global,2);
    error_lower_mat = (regions_year_of_elim_prev_BD_90_mat_lower_median_upper(row_indices_not_global,2) - ...
        regions_year_of_elim_prev_BD_90_mat_lower_median_upper(row_indices_not_global,1));
    error_upper_mat = (regions_year_of_elim_prev_BD_90_mat_lower_median_upper(row_indices_not_global,3) - ...
        regions_year_of_elim_prev_BD_90_mat_lower_median_upper(row_indices_not_global,2));
    errorbar(error_mid_mat,barplot_x_mat,error_lower_mat,error_upper_mat,'r.','horizontal')
    y_mat = regions_year_of_elim_prev_status_quo_mat_lower_median_upper(row_indices_not_global,:);
    dotplots = plot(y_mat(:,2),barplot_x_mat,'o','Color',scenario_colour_mat(status_quo_infant_status_quo_BD_num,:),'MarkerFaceColor',scenario_colour_mat(status_quo_infant_status_quo_BD_num,:));
    for row_num=1:num_regions
        errorbar(y_mat(row_num,2),...
            barplot_x_mat(row_num),...
            y_mat(row_num,2) - y_mat(row_num,1),...
            y_mat(row_num,3) - y_mat(row_num,2),...
            'horizontal','Color',scenario_colour_mat(status_quo_infant_status_quo_BD_num,:))
    end
    xlim([2010 2100])
    line([2030 2030],ylim,'LineStyle',':','Color','k')
    yticklabels(regions_list)
    yticks(1:num_regions)
    xlabel('Year of elimination')
    ylabel('Region')
    lgd = legend([dotplots barplot],ListOfScenarioLabels([status_quo_infant_status_quo_BD_num status_quo_infant_BD_to_90_num]),'Orientation','horizontal','Location','southoutside');
    lgd.NumColumns = 1;
    annotate_subplot(get(gca,'Position'),'a',2)


    clear ListOfScenarioLabels scenario_colour_mat status_quo_infant_status_quo_BD_num status_quo_infant_BD_to_90_num regions_list num_regions row_indices_not_global
    clear regions_year_of_elim_prev_status_quo_mat_lower_median_upper regions_year_of_elim_prev_BD_90_mat_lower_median_upper




    subplot_2 = nexttile;
    regions_list = Fig_4b_data_struc.regions_list;
    num_regions = Fig_4b_data_struc.num_regions;
    region_colour_mat = Fig_4b_data_struc.region_colour_mat;
    regions_BD_coverage_2020_2030_mat = Fig_4b_data_struc.regions_BD_coverage_2020_2030_mat;
    median_regions_year_of_elim_prev_mat = Fig_4b_data_struc.median_regions_year_of_elim_prev_mat;
    lower_regions_year_of_elim_prev_mat = Fig_4b_data_struc.lower_regions_year_of_elim_prev_mat;
    upper_regions_year_of_elim_prev_mat = Fig_4b_data_struc.upper_regions_year_of_elim_prev_mat;


    subplot_2.ColorOrder = region_colour_mat;
    hold on
    region_lines = plot(regions_BD_coverage_2020_2030_mat'/10,NaN(size(median_regions_year_of_elim_prev_mat')),'-','LineWidth',1.5);
    for row_num=1:num_regions
        xvals = regions_BD_coverage_2020_2030_mat(row_num,:)/10;
        yvals = median_regions_year_of_elim_prev_mat(row_num,:);
        colour_val = region_colour_mat(row_num,:);
        index_censored = yvals>2100;
        if sum(index_censored>0)
            assert(sum(index_censored==1))
            assert(yvals(index_censored)==2101) % 2101 indicates a censored value
            assert(yvals(1)==2101)
            plot(xvals(2:end),yvals(2:end),'-','LineWidth',1.5,'Color',colour_val);
        else
            plot(xvals,yvals,'-','LineWidth',1.5,'Color',colour_val);
        end
    end
    for row_num=1:num_regions
        for col_num=1:4
            xval = regions_BD_coverage_2020_2030_mat(row_num,col_num)/10;
            yval_median = median_regions_year_of_elim_prev_mat(row_num,col_num);
            colour_val = region_colour_mat(row_num,:);
            if sum(yval_median>2100)
                continue % do not plot error bar for censored value
            else
                plot(xval,yval_median)
                errorbar(xval,...
                    yval_median,...
                    yval_median - lower_regions_year_of_elim_prev_mat(row_num,col_num),...
                    upper_regions_year_of_elim_prev_mat(row_num,col_num) - yval_median,...
                    'LineWidth',1,'Color',colour_val)
            end
        end
    end
    ylim([2010 2100])
    line(xlim,[2030 2030],'LineStyle',':','Color','k')
    xlabel('Annual increase (% points) in HepB-BD from 2020 to 2030')
    ylabel('Year of elimination')
    lgd = legend(region_lines,regions_list,'Orientation','horizontal','Location','southoutside');
    lgd.NumColumns = 3;
    annotate_subplot(get(gca,'Position'),'b',2)


    clear regions_list num_regions region_colour_mat regions_BD_coverage_2020_2030_mat median_regions_year_of_elim_prev_mat lower_regions_year_of_elim_prev_mat upper_regions_year_of_elim_prev_mat


    assert(strcmp(get(gca, 'FontName'),'Helvetica'))
    subplot_vec = [subplot_1 subplot_2];
    set(subplot_vec,'FontSize',7)


end % end of draw_Figs_4a_4b_fun function





function fig_struct = prepare_Fig_5a_data_struc_fun(...
    num_stochas_runs,...
    cohort_years_vec,cohort_years_pos_vec,...
    scenario_nums_struct,...
    num_regions,regions_list,num_regions_global,regions_global_list,index_global,region_colour_mat,... 
    countries_list,...
    regions_cell_array)


    status_quo_infant_status_quo_BD_num = scenario_nums_struct.status_quo_infant_status_quo_BD_num;
    status_quo_infant_BD_drop_5_num = scenario_nums_struct.status_quo_infant_BD_drop_5_num;
    status_quo_infant_BD_drop_10_num = scenario_nums_struct.status_quo_infant_BD_drop_10_num;
    status_quo_infant_BD_drop_15_num = scenario_nums_struct.status_quo_infant_BD_drop_15_num;
    status_quo_infant_BD_drop_20_num = scenario_nums_struct.status_quo_infant_BD_drop_20_num;


    cohort_years_vec_1yr = 2020:2030;
    num_cohort_years = length(cohort_years_vec_1yr);
    ilower = find(cohort_years_vec>=cohort_years_vec_1yr(1), 1);
    iupper = find(cohort_years_vec>=cohort_years_vec_1yr(end), 1);
    cohort_years_pos_vec_1yr = (cohort_years_pos_vec(ilower)):(cohort_years_pos_vec(iupper));

    [region_excess_deaths_status_quo_BD_drop_5_mat_mean, ...
        region_excess_deaths_status_quo_BD_drop_5_mat_lower_2_5, ...
        region_excess_deaths_status_quo_BD_drop_5_mat_upper_97_5, ...
        region_excess_deaths_status_quo_BD_drop_10_mat_mean, ...
        region_excess_deaths_status_quo_BD_drop_10_mat_lower_2_5, ...
        region_excess_deaths_status_quo_BD_drop_10_mat_upper_97_5, ...
        region_excess_deaths_status_quo_BD_drop_15_mat_mean, ...
        region_excess_deaths_status_quo_BD_drop_15_mat_lower_2_5, ...
        region_excess_deaths_status_quo_BD_drop_15_mat_upper_97_5, ...
        region_excess_deaths_status_quo_BD_drop_20_mat_mean, ...
        region_excess_deaths_status_quo_BD_drop_20_mat_lower_2_5, ...
        region_excess_deaths_status_quo_BD_drop_20_mat_upper_97_5] ...
        = deal(-99 * ones(num_regions_global,1));


    for region_num=1:num_regions_global
    
        region = regions_global_list{region_num};

        [region_excess_deaths_status_quo_BD_drop_5_mat_mean,region_excess_deaths_status_quo_BD_drop_5_mat_lower_2_5,region_excess_deaths_status_quo_BD_drop_5_mat_upper_97_5,...
            ~,~,~,...
            ~,~,~] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_cohort_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD drop 5 2020','Status quo infant & BD',...
            status_quo_infant_BD_drop_5_num,status_quo_infant_status_quo_BD_num,...
            region_excess_deaths_status_quo_BD_drop_5_mat_mean,region_excess_deaths_status_quo_BD_drop_5_mat_lower_2_5,region_excess_deaths_status_quo_BD_drop_5_mat_upper_97_5);

        [region_excess_deaths_status_quo_BD_drop_10_mat_mean,region_excess_deaths_status_quo_BD_drop_10_mat_lower_2_5,region_excess_deaths_status_quo_BD_drop_10_mat_upper_97_5,...
            ~,~,~,...
            ~,~,~] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_cohort_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD drop 10 2020','Status quo infant & BD',...
            status_quo_infant_BD_drop_10_num,status_quo_infant_status_quo_BD_num,...
            region_excess_deaths_status_quo_BD_drop_10_mat_mean,region_excess_deaths_status_quo_BD_drop_10_mat_lower_2_5,region_excess_deaths_status_quo_BD_drop_10_mat_upper_97_5);

        [region_excess_deaths_status_quo_BD_drop_15_mat_mean,region_excess_deaths_status_quo_BD_drop_15_mat_lower_2_5,region_excess_deaths_status_quo_BD_drop_15_mat_upper_97_5,...
            ~,~,~,...
            ~,~,~] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_cohort_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD drop 15 2020','Status quo infant & BD',...
            status_quo_infant_BD_drop_15_num,status_quo_infant_status_quo_BD_num,...
            region_excess_deaths_status_quo_BD_drop_15_mat_mean,region_excess_deaths_status_quo_BD_drop_15_mat_lower_2_5,region_excess_deaths_status_quo_BD_drop_15_mat_upper_97_5);

        [region_excess_deaths_status_quo_BD_drop_20_mat_mean,region_excess_deaths_status_quo_BD_drop_20_mat_lower_2_5,region_excess_deaths_status_quo_BD_drop_20_mat_upper_97_5,...
            ~,~,~,...
            ~,~,~] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_cohort_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD drop 20 2020','Status quo infant & BD',...
            status_quo_infant_BD_drop_20_num,status_quo_infant_status_quo_BD_num,...
            region_excess_deaths_status_quo_BD_drop_20_mat_mean,region_excess_deaths_status_quo_BD_drop_20_mat_lower_2_5,region_excess_deaths_status_quo_BD_drop_20_mat_upper_97_5);

    end % end of region_num for loop


    assert(isequal(size(region_excess_deaths_status_quo_BD_drop_5_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_status_quo_BD_drop_5_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_status_quo_BD_drop_5_mat_upper_97_5),[num_regions_global 1]))
    assert(all(all(region_excess_deaths_status_quo_BD_drop_5_mat_mean>=0)))
    assert(all(all(region_excess_deaths_status_quo_BD_drop_5_mat_lower_2_5>=0)))
    assert(all(all(region_excess_deaths_status_quo_BD_drop_5_mat_upper_97_5>=0)))
    assert(isequal(size(region_excess_deaths_status_quo_BD_drop_10_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_status_quo_BD_drop_10_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_status_quo_BD_drop_10_mat_upper_97_5),[num_regions_global 1]))
    assert(all(all(region_excess_deaths_status_quo_BD_drop_10_mat_mean>=0)))
    assert(all(all(region_excess_deaths_status_quo_BD_drop_10_mat_lower_2_5>=0)))
    assert(all(all(region_excess_deaths_status_quo_BD_drop_10_mat_upper_97_5>=0)))
    assert(isequal(size(region_excess_deaths_status_quo_BD_drop_15_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_status_quo_BD_drop_15_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_status_quo_BD_drop_15_mat_upper_97_5),[num_regions_global 1]))
    assert(all(all(region_excess_deaths_status_quo_BD_drop_15_mat_mean>=0)))
    assert(all(all(region_excess_deaths_status_quo_BD_drop_15_mat_lower_2_5>=0)))
    assert(all(all(region_excess_deaths_status_quo_BD_drop_15_mat_upper_97_5>=0)))
    assert(isequal(size(region_excess_deaths_status_quo_BD_drop_20_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_status_quo_BD_drop_20_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_status_quo_BD_drop_20_mat_upper_97_5),[num_regions_global 1]))
    assert(all(all(region_excess_deaths_status_quo_BD_drop_20_mat_mean>=0)))
    assert(all(all(region_excess_deaths_status_quo_BD_drop_20_mat_lower_2_5>=0)))
    assert(all(all(region_excess_deaths_status_quo_BD_drop_20_mat_upper_97_5>=0)))


    x_vals = categorical({'5% drop','10% drop','15% drop','20% drop'});
    x_vals = reordercats(x_vals,{'5% drop','10% drop','15% drop','20% drop'});

    row_indices_not_global = 1:num_regions_global~=index_global;
    region_excess_deaths_mat_mean = [...
        region_excess_deaths_status_quo_BD_drop_5_mat_mean(row_indices_not_global,:),...
        region_excess_deaths_status_quo_BD_drop_10_mat_mean(row_indices_not_global,:),...
        region_excess_deaths_status_quo_BD_drop_15_mat_mean(row_indices_not_global,:),...
        region_excess_deaths_status_quo_BD_drop_20_mat_mean(row_indices_not_global,:)];
    region_excess_deaths_mat_lower_2_5 = [...
        region_excess_deaths_status_quo_BD_drop_5_mat_lower_2_5(row_indices_not_global,:),...
        region_excess_deaths_status_quo_BD_drop_10_mat_lower_2_5(row_indices_not_global,:),...
        region_excess_deaths_status_quo_BD_drop_15_mat_lower_2_5(row_indices_not_global,:),...
        region_excess_deaths_status_quo_BD_drop_20_mat_lower_2_5(row_indices_not_global,:)];
    region_excess_deaths_mat_upper_97_5 = [...
        region_excess_deaths_status_quo_BD_drop_5_mat_upper_97_5(row_indices_not_global,:),...
        region_excess_deaths_status_quo_BD_drop_10_mat_upper_97_5(row_indices_not_global,:),...
        region_excess_deaths_status_quo_BD_drop_15_mat_upper_97_5(row_indices_not_global,:),...
        region_excess_deaths_status_quo_BD_drop_20_mat_upper_97_5(row_indices_not_global,:)];
    assert(isequal(size(region_excess_deaths_mat_mean),[num_regions 4]))
    assert(isequal(size(region_excess_deaths_mat_lower_2_5),[num_regions 4]))
    assert(isequal(size(region_excess_deaths_mat_upper_97_5),[num_regions 4]))


    fig_struct.regions_list = regions_list;
    fig_struct.num_regions = num_regions;
    fig_struct.region_colour_mat = region_colour_mat;
    fig_struct.x_vals = x_vals;
    fig_struct.region_excess_deaths_mat_mean = region_excess_deaths_mat_mean;
    fig_struct.region_excess_deaths_mat_lower_2_5 = region_excess_deaths_mat_lower_2_5;
    fig_struct.region_excess_deaths_mat_upper_97_5 = region_excess_deaths_mat_upper_97_5;


    disp('Additional HBV-related deaths globally in scenario 3a relative to scenario 1 in the 2020 and 2030 birth cohorts')
    mean_num = region_excess_deaths_status_quo_BD_drop_5_mat_mean(index_global,:);
    lower_num = region_excess_deaths_status_quo_BD_drop_5_mat_lower_2_5(index_global,:);
    upper_num = region_excess_deaths_status_quo_BD_drop_5_mat_upper_97_5(index_global,:);
    assert(mean_num>lower_num)
    assert(mean_num<upper_num)
    disp([num2str(round(mean_num)) ' (' num2str(round(lower_num)) ', ' num2str(round(upper_num)) ')'])

    disp('Additional HBV-related deaths globally in scenario 3b relative to scenario 1 in the 2020 and 2030 birth cohorts')
    mean_num = region_excess_deaths_status_quo_BD_drop_10_mat_mean(index_global,:);
    lower_num = region_excess_deaths_status_quo_BD_drop_10_mat_lower_2_5(index_global,:);
    upper_num = region_excess_deaths_status_quo_BD_drop_10_mat_upper_97_5(index_global,:);
    assert(mean_num>lower_num)
    assert(mean_num<upper_num)
    disp([num2str(round(mean_num)) ' (' num2str(round(lower_num)) ', ' num2str(round(upper_num)) ')'])

    disp('Additional HBV-related deaths globally in scenario 3c relative to scenario 1 in the 2020 and 2030 birth cohorts')
    mean_num = region_excess_deaths_status_quo_BD_drop_15_mat_mean(index_global,:);
    lower_num = region_excess_deaths_status_quo_BD_drop_15_mat_lower_2_5(index_global,:);
    upper_num = region_excess_deaths_status_quo_BD_drop_15_mat_upper_97_5(index_global,:);
    assert(mean_num>lower_num)
    assert(mean_num<upper_num)
    disp([num2str(round(mean_num)) ' (' num2str(round(lower_num)) ', ' num2str(round(upper_num)) ')'])

    disp('Additional HBV-related deaths globally in scenario 3d relative to scenario 1 in the 2020 and 2030 birth cohorts')
    mean_num = region_excess_deaths_status_quo_BD_drop_20_mat_mean(index_global,:);
    lower_num = region_excess_deaths_status_quo_BD_drop_20_mat_lower_2_5(index_global,:);
    upper_num = region_excess_deaths_status_quo_BD_drop_20_mat_upper_97_5(index_global,:);
    assert(mean_num>lower_num)
    assert(mean_num<upper_num)
    disp([num2str(round(mean_num)) ' (' num2str(round(lower_num)) ', ' num2str(round(upper_num)) ')'])


end % end of prepare_Fig_5a_data_struc_fun function






function fig_struct = prepare_Fig_5b_data_struc_fun(...
    num_stochas_runs,...
    cohort_years_vec,cohort_years_pos_vec,...
    scenario_nums_struct,...
    num_regions,regions_list,num_regions_global,regions_global_list,index_global,region_colour_mat,...
    countries_list,...
    regions_cell_array)


    status_quo_infant_BD_to_90_num = scenario_nums_struct.status_quo_infant_BD_to_90_num;
    BD_delayed_2023_2030_num = scenario_nums_struct.BD_delayed_2023_2030_num;
    BD_delayed_slowed_2023_2033_num = scenario_nums_struct.BD_delayed_slowed_2023_2033_num;
    BD_delayed_slowed_2025_2040_num = scenario_nums_struct.BD_delayed_slowed_2025_2040_num;


    cohort_years_vec_1yr = 2020:2030;
    num_cohort_years = length(cohort_years_vec_1yr);
    ilower = find(cohort_years_vec>=cohort_years_vec_1yr(1), 1);
    iupper = find(cohort_years_vec>=cohort_years_vec_1yr(end), 1);
    cohort_years_pos_vec_1yr = (cohort_years_pos_vec(ilower)):(cohort_years_pos_vec(iupper));

    [region_excess_deaths_BD_90_2023_2030_mat_mean, ...
        region_excess_deaths_BD_90_2023_2030_mat_lower_2_5, ...
        region_excess_deaths_BD_90_2023_2030_mat_upper_97_5, ...
        region_excess_deaths_BD_90_2023_2033_mat_mean, ...
        region_excess_deaths_BD_90_2023_2033_mat_lower_2_5, ...
        region_excess_deaths_BD_90_2023_2033_mat_upper_97_5, ...
        region_excess_deaths_BD_90_2025_2040_mat_mean, ...
        region_excess_deaths_BD_90_2025_2040_mat_lower_2_5, ...
        region_excess_deaths_BD_90_2025_2040_mat_upper_97_5] ...
        = deal(-99 * ones(num_regions_global,1));


    for region_num=1:num_regions_global
        
        region = regions_global_list{region_num};

        [region_excess_deaths_BD_90_2023_2030_mat_mean,region_excess_deaths_BD_90_2023_2030_mat_lower_2_5,region_excess_deaths_BD_90_2023_2030_mat_upper_97_5,...
            ~,~,~,...
            ~,~,~] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_cohort_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD delayed expansion 2023 to 2030','Status quo infant & BD expansion to 90%',...
            BD_delayed_2023_2030_num,status_quo_infant_BD_to_90_num,...
            region_excess_deaths_BD_90_2023_2030_mat_mean,region_excess_deaths_BD_90_2023_2030_mat_lower_2_5,region_excess_deaths_BD_90_2023_2030_mat_upper_97_5);

        [region_excess_deaths_BD_90_2023_2033_mat_mean,region_excess_deaths_BD_90_2023_2033_mat_lower_2_5,region_excess_deaths_BD_90_2023_2033_mat_upper_97_5,...
            ~,~,~,...
            ~,~,~] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_cohort_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD delayed expansion 2023 to 2033','Status quo infant & BD expansion to 90%',...
            BD_delayed_slowed_2023_2033_num,status_quo_infant_BD_to_90_num,...
            region_excess_deaths_BD_90_2023_2033_mat_mean,region_excess_deaths_BD_90_2023_2033_mat_lower_2_5,region_excess_deaths_BD_90_2023_2033_mat_upper_97_5);

        [region_excess_deaths_BD_90_2025_2040_mat_mean,region_excess_deaths_BD_90_2025_2040_mat_lower_2_5,region_excess_deaths_BD_90_2025_2040_mat_upper_97_5,...
            ~,~,~,...
            ~,~,~] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_cohort_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD delayed expansion 2025 to 2040','Status quo infant & BD expansion to 90%',...
            BD_delayed_slowed_2025_2040_num,status_quo_infant_BD_to_90_num,...
            region_excess_deaths_BD_90_2025_2040_mat_mean,region_excess_deaths_BD_90_2025_2040_mat_lower_2_5,region_excess_deaths_BD_90_2025_2040_mat_upper_97_5);

    end % end of region_num for loop


    assert(isequal(size(region_excess_deaths_BD_90_2023_2030_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_BD_90_2023_2030_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_BD_90_2023_2030_mat_upper_97_5),[num_regions_global 1]))
    assert(all(all(region_excess_deaths_BD_90_2023_2030_mat_mean>=0)))
    assert(all(all(region_excess_deaths_BD_90_2023_2030_mat_lower_2_5>=0)))
    assert(all(all(region_excess_deaths_BD_90_2023_2030_mat_upper_97_5>=0)))
    assert(isequal(size(region_excess_deaths_BD_90_2023_2033_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_BD_90_2023_2033_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_BD_90_2023_2033_mat_upper_97_5),[num_regions_global 1]))
    assert(all(all(region_excess_deaths_BD_90_2023_2033_mat_mean>=0)))
    assert(all(all(region_excess_deaths_BD_90_2023_2033_mat_lower_2_5>=0)))
    assert(all(all(region_excess_deaths_BD_90_2023_2033_mat_upper_97_5>=0)))
    assert(isequal(size(region_excess_deaths_BD_90_2025_2040_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_BD_90_2025_2040_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_BD_90_2025_2040_mat_upper_97_5),[num_regions_global 1]))
    assert(all(all(region_excess_deaths_BD_90_2025_2040_mat_mean>=0)))
    assert(all(all(region_excess_deaths_BD_90_2025_2040_mat_lower_2_5>=0)))
    assert(all(all(region_excess_deaths_BD_90_2025_2040_mat_upper_97_5>=0)))


    x_vals = categorical({'2023 to 2030','2023 to 2033','2025 to 2040'});
    x_vals = reordercats(x_vals,{'2023 to 2030','2023 to 2033','2025 to 2040'});

    row_indices_not_global = 1:num_regions_global~=index_global;
    region_excess_deaths_mat_mean = [...
        region_excess_deaths_BD_90_2023_2030_mat_mean(row_indices_not_global,:),...
        region_excess_deaths_BD_90_2023_2033_mat_mean(row_indices_not_global,:),...
        region_excess_deaths_BD_90_2025_2040_mat_mean(row_indices_not_global,:)];
    region_excess_deaths_mat_lower_2_5 = [...
        region_excess_deaths_BD_90_2023_2030_mat_lower_2_5(row_indices_not_global,:),...
        region_excess_deaths_BD_90_2023_2033_mat_lower_2_5(row_indices_not_global,:),...
        region_excess_deaths_BD_90_2025_2040_mat_lower_2_5(row_indices_not_global,:)];
    region_excess_deaths_mat_upper_97_5 = [...
        region_excess_deaths_BD_90_2023_2030_mat_upper_97_5(row_indices_not_global,:),...
        region_excess_deaths_BD_90_2023_2033_mat_upper_97_5(row_indices_not_global,:),...
        region_excess_deaths_BD_90_2025_2040_mat_upper_97_5(row_indices_not_global,:)];
    assert(isequal(size(region_excess_deaths_mat_mean),[num_regions 3]))
    assert(isequal(size(region_excess_deaths_mat_lower_2_5),[num_regions 3]))
    assert(isequal(size(region_excess_deaths_mat_upper_97_5),[num_regions 3]))


    fig_struct.regions_list = regions_list;
    fig_struct.num_regions = num_regions;
    fig_struct.region_colour_mat = region_colour_mat;
    fig_struct.x_vals = x_vals;
    fig_struct.region_excess_deaths_mat_mean = region_excess_deaths_mat_mean;
    fig_struct.region_excess_deaths_mat_lower_2_5 = region_excess_deaths_mat_lower_2_5;
    fig_struct.region_excess_deaths_mat_upper_97_5 = region_excess_deaths_mat_upper_97_5;


    disp('Additional HBV-related deaths globally in scenario 4a relative to scenario 2 in the 2020 and 2030 birth cohorts')
    mean_num = region_excess_deaths_BD_90_2023_2030_mat_mean(index_global,:);
    lower_num = region_excess_deaths_BD_90_2023_2030_mat_lower_2_5(index_global,:);
    upper_num = region_excess_deaths_BD_90_2023_2030_mat_upper_97_5(index_global,:);
    assert(mean_num>lower_num)
    assert(mean_num<upper_num)
    disp([num2str(round(mean_num)) ' (' num2str(round(lower_num)) ', ' num2str(round(upper_num)) ')'])

    disp('Additional HBV-related deaths globally in scenario 4b relative to scenario 2 in the 2020 and 2030 birth cohorts')
    mean_num = region_excess_deaths_BD_90_2023_2033_mat_mean(index_global,:);
    lower_num = region_excess_deaths_BD_90_2023_2033_mat_lower_2_5(index_global,:);
    upper_num = region_excess_deaths_BD_90_2023_2033_mat_upper_97_5(index_global,:);
    assert(mean_num>lower_num)
    assert(mean_num<upper_num)
    disp([num2str(round(mean_num)) ' (' num2str(round(lower_num)) ', ' num2str(round(upper_num)) ')'])

    disp('Additional HBV-related deaths globally in scenario 4c relative to scenario 2 in the 2020 and 2030 birth cohorts')
    mean_num = region_excess_deaths_BD_90_2025_2040_mat_mean(index_global,:);
    lower_num = region_excess_deaths_BD_90_2025_2040_mat_lower_2_5(index_global,:);
    upper_num = region_excess_deaths_BD_90_2025_2040_mat_upper_97_5(index_global,:);
    assert(mean_num>lower_num)
    assert(mean_num<upper_num)
    disp([num2str(round(mean_num)) ' (' num2str(round(lower_num)) ', ' num2str(round(upper_num)) ')'])


end % end of prepare_Fig_5b_data_struc_fun function





function draw_Figs_5a_5b_fun(Fig_5a_data_struc,Fig_5b_data_struc)


    f = figure('Units','centimeters');
    f.Position(3:4) = [17.5 8];
    tiledlayout(1, 2)


    subplot_1 = nexttile;
    regions_list = Fig_5a_data_struc.regions_list;
    num_regions = Fig_5a_data_struc.num_regions;
    region_colour_mat = Fig_5a_data_struc.region_colour_mat;
    x_vals = Fig_5a_data_struc.x_vals;
    region_excess_deaths_mat_mean = Fig_5a_data_struc.region_excess_deaths_mat_mean;
    region_excess_deaths_mat_lower_2_5 = Fig_5a_data_struc.region_excess_deaths_mat_lower_2_5;
    region_excess_deaths_mat_upper_97_5 = Fig_5a_data_struc.region_excess_deaths_mat_upper_97_5;


    the_denominator = 1e3;
    ylabel_str = 'thousands';

    subplot_1.ColorOrder = region_colour_mat;
    hold on
    barplot = bar(x_vals,region_excess_deaths_mat_mean'/the_denominator,'stacked');
    error_mid_mat = cumsum(region_excess_deaths_mat_mean',2)/the_denominator;
    error_lower_mat = (region_excess_deaths_mat_mean - region_excess_deaths_mat_lower_2_5)'/the_denominator;
    error_upper_mat = (region_excess_deaths_mat_upper_97_5 - region_excess_deaths_mat_mean)'/the_denominator;
    for region_num=1:num_regions
        errorbar(x_vals,error_mid_mat(:,region_num),error_lower_mat(:,region_num),error_upper_mat(:,region_num),'.k')
    end
    ylim([0 22e3]/the_denominator)
    xlabel('Scenario')
    ylabel(['Additional HBV-related deaths (' ylabel_str ')'])
    lgd = legend(regions_list,'Location','southoutside','Orientation','horizontal');
    lgd.NumColumns = 3;
    annotate_subplot(get(gca,'Position'),'a',2)


    clear regions_list num_regions region_colour_mat x_vals region_excess_deaths_mat_mean region_excess_deaths_mat_lower_2_5 region_excess_deaths_mat_upper_97_5




    subplot_2 = nexttile;
    regions_list = Fig_5b_data_struc.regions_list;
    num_regions = Fig_5b_data_struc.num_regions;
    region_colour_mat = Fig_5b_data_struc.region_colour_mat;
    x_vals = Fig_5b_data_struc.x_vals;
    region_excess_deaths_mat_mean = Fig_5b_data_struc.region_excess_deaths_mat_mean;
    region_excess_deaths_mat_lower_2_5 = Fig_5b_data_struc.region_excess_deaths_mat_lower_2_5;
    region_excess_deaths_mat_upper_97_5 = Fig_5b_data_struc.region_excess_deaths_mat_upper_97_5;


    the_denominator = 1e3;
    ylabel_str = 'thousands';

    subplot_2.ColorOrder = region_colour_mat;
    hold on
    barplot = bar(x_vals,region_excess_deaths_mat_mean'/the_denominator,'stacked');
    error_mid_mat = cumsum(region_excess_deaths_mat_mean',2)/the_denominator;
    error_lower_mat = (region_excess_deaths_mat_mean - region_excess_deaths_mat_lower_2_5)'/the_denominator;
    error_upper_mat = (region_excess_deaths_mat_upper_97_5 - region_excess_deaths_mat_mean)'/the_denominator;
    for region_num=1:num_regions
        errorbar(x_vals,error_mid_mat(:,region_num),error_lower_mat(:,region_num),error_upper_mat(:,region_num),'.k')
    end
    ylim([0 700e3]/the_denominator)
    xlabel('Scenario')
    ylabel(['Additional HBV-related deaths (' ylabel_str ')'])
    lgd = legend(regions_list,'Location','southoutside','Orientation','horizontal');
    lgd.NumColumns = 3;
    annotate_subplot(get(gca,'Position'),'b',2)


    clear regions_list num_regions region_colour_mat x_vals region_excess_deaths_mat_mean region_excess_deaths_mat_lower_2_5 region_excess_deaths_mat_upper_97_5


    assert(strcmp(get(gca, 'FontName'),'Helvetica'))
    subplot_vec = [subplot_1 subplot_2];
    set(subplot_vec,'FontSize',7)


end % end of draw_Figs_5a_5b_fun function






function draw_Suppl_Tables_3_4_fun(...
    num_stochas_runs,...
    calandar_years_vec,calandar_years_pos_vec,...
    cohort_years_vec,cohort_years_pos_vec,...
    scenario_nums_struct,...
    num_regions_global,regions_global_list,...
    countries_list,...
    regions_cell_array)


    status_quo_infant_status_quo_BD_num = scenario_nums_struct.status_quo_infant_status_quo_BD_num;
    status_quo_infant_BD_to_25_num = scenario_nums_struct.status_quo_infant_BD_to_25_num;
    status_quo_infant_BD_to_90_num = scenario_nums_struct.status_quo_infant_BD_to_90_num;
    status_quo_infant_BD_drop_5_num = scenario_nums_struct.status_quo_infant_BD_drop_5_num;
    status_quo_infant_BD_drop_20_num = scenario_nums_struct.status_quo_infant_BD_drop_20_num;
    BD_delayed_2023_2030_num = scenario_nums_struct.BD_delayed_2023_2030_num;
    BD_delayed_slowed_2025_2040_num = scenario_nums_struct.BD_delayed_slowed_2025_2040_num;


    cohort_years_vec_1yr = 2020:2030;
    num_years = length(cohort_years_vec_1yr);
    ilower = find(cohort_years_vec>=cohort_years_vec_1yr(1), 1);
    iupper = find(cohort_years_vec>=cohort_years_vec_1yr(end), 1);
    cohort_years_pos_vec_1yr = (cohort_years_pos_vec(ilower)):(cohort_years_pos_vec(iupper));


    [regions_year_of_elim_prev_1_mat_lower_median_upper,...
        regions_year_of_elim_prev_2a_mat_lower_median_upper,...
        regions_year_of_elim_prev_2d_mat_lower_median_upper,...
        regions_year_of_elim_prev_3a_mat_lower_median_upper,...
        regions_year_of_elim_prev_3d_mat_lower_median_upper,...
        regions_year_of_elim_prev_4a_mat_lower_median_upper,...
        regions_year_of_elim_prev_4c_mat_lower_median_upper] ...
        = deal(zeros(num_regions_global,3));


    [region_excess_deaths_1_2a_mat_mean, ...
        region_excess_deaths_1_2a_mat_lower_2_5, ...
        region_excess_deaths_1_2a_mat_upper_97_5, ...
        region_excess_deaths_1_2d_mat_mean, ...
        region_excess_deaths_1_2d_mat_lower_2_5, ...
        region_excess_deaths_1_2d_mat_upper_97_5, ...
        region_excess_deaths_1_3a_mat_mean, ...
        region_excess_deaths_1_3a_mat_lower_2_5, ...
        region_excess_deaths_1_3a_mat_upper_97_5, ...
        region_excess_deaths_1_3d_mat_mean, ...
        region_excess_deaths_1_3d_mat_lower_2_5, ...
        region_excess_deaths_1_3d_mat_upper_97_5, ...
        region_excess_deaths_1_4a_mat_mean, ...
        region_excess_deaths_1_4a_mat_lower_2_5, ...
        region_excess_deaths_1_4a_mat_upper_97_5, ...
        region_excess_deaths_1_4c_mat_mean, ...
        region_excess_deaths_1_4c_mat_lower_2_5, ...
        region_excess_deaths_1_4c_mat_upper_97_5] ... 
        = deal(-99 * ones(num_regions_global,1));


    for region_num=1:num_regions_global
        
        region = regions_global_list{region_num};

        regions_year_of_elim_prev_1_mat_lower_median_upper = ...
            find_year_of_elim_fun(countries_list,calandar_years_vec,calandar_years_pos_vec,...
            'regions',region,region_num,regions_cell_array,num_stochas_runs,0.1,...
            'Status quo infant & BD',status_quo_infant_status_quo_BD_num,...
            regions_year_of_elim_prev_1_mat_lower_median_upper);

        regions_year_of_elim_prev_2a_mat_lower_median_upper = ...
            find_year_of_elim_fun(countries_list,calandar_years_vec,calandar_years_pos_vec,...
            'regions',region,region_num,regions_cell_array,num_stochas_runs,0.1,...
            'Status quo infant & BD expansion to 25%',status_quo_infant_BD_to_25_num,...
            regions_year_of_elim_prev_2a_mat_lower_median_upper);

        regions_year_of_elim_prev_2d_mat_lower_median_upper = ...
            find_year_of_elim_fun(countries_list,calandar_years_vec,calandar_years_pos_vec,...
            'regions',region,region_num,regions_cell_array,num_stochas_runs,0.1,...
            'Status quo infant & BD expansion to 90%',status_quo_infant_BD_to_90_num,...
            regions_year_of_elim_prev_2d_mat_lower_median_upper);

        regions_year_of_elim_prev_3a_mat_lower_median_upper = ...
            find_year_of_elim_fun(countries_list,calandar_years_vec,calandar_years_pos_vec,...
            'regions',region,region_num,regions_cell_array,num_stochas_runs,0.1,...
            'Status quo infant & BD drop 5 2020',status_quo_infant_BD_drop_5_num,...
            regions_year_of_elim_prev_3a_mat_lower_median_upper);

        regions_year_of_elim_prev_3d_mat_lower_median_upper = ...
            find_year_of_elim_fun(countries_list,calandar_years_vec,calandar_years_pos_vec,...
            'regions',region,region_num,regions_cell_array,num_stochas_runs,0.1,...
            'Status quo infant & BD drop 20 2020',status_quo_infant_BD_drop_20_num,...
            regions_year_of_elim_prev_3d_mat_lower_median_upper);

        regions_year_of_elim_prev_4a_mat_lower_median_upper = ...
            find_year_of_elim_fun(countries_list,calandar_years_vec,calandar_years_pos_vec,...
            'regions',region,region_num,regions_cell_array,num_stochas_runs,0.1,...
            'Status quo infant & BD delayed expansion 2023 to 2030',BD_delayed_2023_2030_num,...
            regions_year_of_elim_prev_4a_mat_lower_median_upper);

        regions_year_of_elim_prev_4c_mat_lower_median_upper = ...
            find_year_of_elim_fun(countries_list,calandar_years_vec,calandar_years_pos_vec,...
            'regions',region,region_num,regions_cell_array,num_stochas_runs,0.1,...
            'Status quo infant & BD delayed expansion 2025 to 2040',BD_delayed_slowed_2025_2040_num,...
            regions_year_of_elim_prev_4c_mat_lower_median_upper);

        [region_excess_deaths_1_2a_mat_mean,region_excess_deaths_1_2a_mat_lower_2_5,region_excess_deaths_1_2a_mat_upper_97_5] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD','Status quo infant & BD expansion to 25%',...
            status_quo_infant_status_quo_BD_num,status_quo_infant_BD_to_25_num,...
            region_excess_deaths_1_2a_mat_mean,region_excess_deaths_1_2a_mat_lower_2_5,region_excess_deaths_1_2a_mat_upper_97_5);

        [region_excess_deaths_1_2d_mat_mean,region_excess_deaths_1_2d_mat_lower_2_5,region_excess_deaths_1_2d_mat_upper_97_5] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD','Status quo infant & BD expansion to 90%',...
            status_quo_infant_status_quo_BD_num,status_quo_infant_BD_to_90_num,...
            region_excess_deaths_1_2d_mat_mean,region_excess_deaths_1_2d_mat_lower_2_5,region_excess_deaths_1_2d_mat_upper_97_5);

        [region_excess_deaths_1_3a_mat_mean,region_excess_deaths_1_3a_mat_lower_2_5,region_excess_deaths_1_3a_mat_upper_97_5] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD','Status quo infant & BD drop 5 2020',...
            status_quo_infant_status_quo_BD_num,status_quo_infant_BD_drop_5_num,...
            region_excess_deaths_1_3a_mat_mean,region_excess_deaths_1_3a_mat_lower_2_5,region_excess_deaths_1_3a_mat_upper_97_5);

        [region_excess_deaths_1_3d_mat_mean,region_excess_deaths_1_3d_mat_lower_2_5,region_excess_deaths_1_3d_mat_upper_97_5] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD','Status quo infant & BD drop 20 2020',...
            status_quo_infant_status_quo_BD_num,status_quo_infant_BD_drop_20_num,...
            region_excess_deaths_1_3d_mat_mean,region_excess_deaths_1_3d_mat_lower_2_5,region_excess_deaths_1_3d_mat_upper_97_5);

        [region_excess_deaths_1_4a_mat_mean,region_excess_deaths_1_4a_mat_lower_2_5,region_excess_deaths_1_4a_mat_upper_97_5] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD','Status quo infant & BD delayed expansion 2023 to 2030',...
            status_quo_infant_status_quo_BD_num,BD_delayed_2023_2030_num,...
            region_excess_deaths_1_4a_mat_mean,region_excess_deaths_1_4a_mat_lower_2_5,region_excess_deaths_1_4a_mat_upper_97_5);

        [region_excess_deaths_1_4c_mat_mean,region_excess_deaths_1_4c_mat_lower_2_5,region_excess_deaths_1_4c_mat_upper_97_5] = ...
            difference_in_deaths_fun(true,...
            countries_list,cohort_years_vec_1yr,num_years,cohort_years_pos_vec_1yr,'regions',region,region_num,regions_cell_array,num_stochas_runs,'Incid_Deaths_1yr_approx_birth_cohorts',...
            'Status quo infant & BD','Status quo infant & BD delayed expansion 2025 to 2040',...
            status_quo_infant_status_quo_BD_num,BD_delayed_slowed_2025_2040_num,...
            region_excess_deaths_1_4c_mat_mean,region_excess_deaths_1_4c_mat_lower_2_5,region_excess_deaths_1_4c_mat_upper_97_5);

    end % end of region_num for loop


    assert(isequal(size(regions_year_of_elim_prev_1_mat_lower_median_upper),[num_regions_global 3]))
    assert(isequal(size(regions_year_of_elim_prev_2a_mat_lower_median_upper),[num_regions_global 3]))
    assert(isequal(size(regions_year_of_elim_prev_2d_mat_lower_median_upper),[num_regions_global 3]))
    assert(isequal(size(regions_year_of_elim_prev_3a_mat_lower_median_upper),[num_regions_global 3]))
    assert(isequal(size(regions_year_of_elim_prev_3d_mat_lower_median_upper),[num_regions_global 3]))
    assert(isequal(size(regions_year_of_elim_prev_4a_mat_lower_median_upper),[num_regions_global 3]))
    assert(isequal(size(regions_year_of_elim_prev_4c_mat_lower_median_upper),[num_regions_global 3]))
    assert(all(all(regions_year_of_elim_prev_1_mat_lower_median_upper>=0)))
    assert(all(all(regions_year_of_elim_prev_2a_mat_lower_median_upper>=0)))
    assert(all(all(regions_year_of_elim_prev_2d_mat_lower_median_upper>=0)))
    assert(all(all(regions_year_of_elim_prev_3a_mat_lower_median_upper>=0)))
    assert(all(all(regions_year_of_elim_prev_3d_mat_lower_median_upper>=0)))
    assert(all(all(regions_year_of_elim_prev_4a_mat_lower_median_upper>=0)))
    assert(all(all(regions_year_of_elim_prev_4c_mat_lower_median_upper>=0)))
    assert(isequal(size(region_excess_deaths_1_2a_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_2a_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_2a_mat_upper_97_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_2d_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_2d_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_2d_mat_upper_97_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_3a_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_3a_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_3a_mat_upper_97_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_3d_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_3d_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_3d_mat_upper_97_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_4a_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_4a_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_4a_mat_upper_97_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_4c_mat_mean),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_4c_mat_lower_2_5),[num_regions_global 1]))
    assert(isequal(size(region_excess_deaths_1_4c_mat_upper_97_5),[num_regions_global 1]))
    assert(all(all(region_excess_deaths_1_2a_mat_mean>=0)))
    assert(all(all(region_excess_deaths_1_2a_mat_lower_2_5>=0)))
    assert(all(all(region_excess_deaths_1_2a_mat_upper_97_5>=0)))
    assert(all(all(region_excess_deaths_1_2d_mat_mean>=0)))
    assert(all(all(region_excess_deaths_1_2d_mat_lower_2_5>=0)))
    assert(all(all(region_excess_deaths_1_2d_mat_upper_97_5>=0)))
    assert(all(all(region_excess_deaths_1_3a_mat_mean<=0)))
    assert(all(all(region_excess_deaths_1_3a_mat_lower_2_5<=0)))
    assert(all(all(region_excess_deaths_1_3a_mat_upper_97_5<=0)))
    assert(all(all(region_excess_deaths_1_3d_mat_mean<=0)))
    assert(all(all(region_excess_deaths_1_3d_mat_lower_2_5<=0)))
    assert(all(all(region_excess_deaths_1_3d_mat_upper_97_5<=0)))
    assert(all(all(region_excess_deaths_1_4a_mat_mean>=0)))
    assert(all(all(region_excess_deaths_1_4a_mat_lower_2_5>=0)))
    assert(all(all(region_excess_deaths_1_4a_mat_upper_97_5>=0)))
    assert(all(all(region_excess_deaths_1_4c_mat_mean>=0)))
    assert(all(all(region_excess_deaths_1_4c_mat_lower_2_5>=0)))
    assert(all(all(region_excess_deaths_1_4c_mat_upper_97_5>=0)))



    filename = 'Suppl_tables_3_4.xlsx';


    tab_mean=table(...
        regions_global_list',...
        round(regions_year_of_elim_prev_1_mat_lower_median_upper(:,2)), ...
        round(regions_year_of_elim_prev_2a_mat_lower_median_upper(:,2)), ...
        round(regions_year_of_elim_prev_2d_mat_lower_median_upper(:,2)), ...
        round(regions_year_of_elim_prev_3a_mat_lower_median_upper(:,2)), ...
        round(regions_year_of_elim_prev_3d_mat_lower_median_upper(:,2)), ...
        round(regions_year_of_elim_prev_4a_mat_lower_median_upper(:,2)), ...
        round(regions_year_of_elim_prev_4c_mat_lower_median_upper(:,2)) ...
        );


    var_names = {...
        'Region', ...
        'Year_of_elim_1_median', ...
        'Year_of_elim_2a_median', ...
        'Year_of_elim_2d_median', ...
        'Year_of_elim_3a_median', ...
        'Year_of_elim_3d_median', ...
        'Year_of_elim_4a_median', ...
        'Year_of_elim_4c_median' ...
        };
    num_varnames = length(var_names);
    tab_mean.Properties.VariableNames = var_names;


    for varname_num = 2:num_varnames
        curr_varname = var_names{varname_num};
        temp_var = tab_mean.(curr_varname);
        temp_var = string(temp_var);
        temp_var(strcmp(temp_var,'1980')) = {'before 1980'};
        temp_var(strcmp(temp_var,'2101')) = {'after 2100'};
        tab_mean.(curr_varname) = temp_var;
    end


    %writetable(tab_mean,filename,'Sheet','Region Year of elimination median')


    tab_lower=table(...
        regions_global_list',...
        round(regions_year_of_elim_prev_1_mat_lower_median_upper(:,1)), ...
        round(regions_year_of_elim_prev_2a_mat_lower_median_upper(:,1)), ...
        round(regions_year_of_elim_prev_2d_mat_lower_median_upper(:,1)), ...
        round(regions_year_of_elim_prev_3a_mat_lower_median_upper(:,1)), ...
        round(regions_year_of_elim_prev_3d_mat_lower_median_upper(:,1)), ...
        round(regions_year_of_elim_prev_4a_mat_lower_median_upper(:,1)), ...
        round(regions_year_of_elim_prev_4c_mat_lower_median_upper(:,1)) ...
        );


    var_names = {...
        'Region', ...
        'Year_of_elim_1_lower_2_5', ...
        'Year_of_elim_2a_lower_2_5', ...
        'Year_of_elim_2d_lower_2_5', ...
        'Year_of_elim_3a_lower_2_5', ...
        'Year_of_elim_3d_lower_2_5', ...
        'Year_of_elim_4a_lower_2_5', ...
        'Year_of_elim_4c_lower_2_5' ...
        };
    num_varnames = length(var_names);
    tab_lower.Properties.VariableNames = var_names;


    for varname_num = 2:num_varnames
        curr_varname = var_names{varname_num};
        temp_var = tab_lower.(curr_varname);
        temp_var = string(temp_var);
        temp_var(strcmp(temp_var,'1980')) = {'before 1980'};
        temp_var(strcmp(temp_var,'2101')) = {'after 2100'};
        tab_lower.(curr_varname) = temp_var;
    end


    %writetable(tab_lower,filename,'Sheet','Region Year of elimination lower 2.5')


    tab_upper=table(...
        regions_global_list',...
        round(regions_year_of_elim_prev_1_mat_lower_median_upper(:,3)), ...
        round(regions_year_of_elim_prev_2a_mat_lower_median_upper(:,3)), ...
        round(regions_year_of_elim_prev_2d_mat_lower_median_upper(:,3)), ...
        round(regions_year_of_elim_prev_3a_mat_lower_median_upper(:,3)), ...
        round(regions_year_of_elim_prev_3d_mat_lower_median_upper(:,3)), ...
        round(regions_year_of_elim_prev_4a_mat_lower_median_upper(:,3)), ...
        round(regions_year_of_elim_prev_4c_mat_lower_median_upper(:,3)) ...
        );


    var_names = {...
        'Region', ...
        'Year_of_elim_1_upper_97_5', ...
        'Year_of_elim_2a_upper_97_5', ...
        'Year_of_elim_2d_upper_97_5', ...
        'Year_of_elim_3a_upper_97_5', ...
        'Year_of_elim_3d_upper_97_5', ...
        'Year_of_elim_4a_upper_97_5', ...
        'Year_of_elim_4c_upper_97_5' ...
        };
    num_varnames = length(var_names);
    tab_upper.Properties.VariableNames = var_names;


    for varname_num = 2:num_varnames
        curr_varname = var_names{varname_num};
        temp_var = tab_upper.(curr_varname);
        temp_var = string(temp_var);
        temp_var(strcmp(temp_var,'1980')) = {'before 1980'};
        temp_var(strcmp(temp_var,'2101')) = {'after 2100'};
        tab_upper.(curr_varname) = temp_var;
    end


    %writetable(tab_upper,filename,'Sheet','Region Year of elimination upper 97.5')


    tab = build_mean_min_max_string(tab_mean,tab_lower,tab_upper,num_regions_global,8);


    var_names = {...
        'Region', ...
        'Year_of_elim_1', ...
        'Year_of_elim_2a', ...
        'Year_of_elim_2d', ...
        'Year_of_elim_3a', ...
        'Year_of_elim_3d', ...
        'Year_of_elim_4a', ...
        'Year_of_elim_4c' ...
        };
    num_varnames = length(var_names);
    tab.Properties.VariableNames = var_names;


    writetable(tab,filename,'Sheet','Region Year of elimination')


    empty_val_cell = repmat({'-'},num_regions_global,1);

    tab_mean=table(...
        regions_global_list',...
        empty_val_cell, ...
        round(region_excess_deaths_1_2a_mat_mean), ...
        round(region_excess_deaths_1_2d_mat_mean), ...
        round(region_excess_deaths_1_3a_mat_mean), ...
        round(region_excess_deaths_1_3d_mat_mean), ...
        round(region_excess_deaths_1_4a_mat_mean), ...
        round(region_excess_deaths_1_4c_mat_mean) ...
        );


    var_names = {...
        'Region', ...
        'Excess_deaths_1_mean', ...
        'Excess_deaths_2a_mean', ...
        'Excess_deaths_2d_mean', ...
        'Excess_deaths_3a_mean', ...
        'Excess_deaths_3d_mean', ...
        'Excess_deaths_4a_mean', ...
        'Excess_deaths_4c_mean' ...
        };
    num_varnames = length(var_names);
    tab_mean.Properties.VariableNames = var_names;


    for varname_num = 3:num_varnames
        curr_varname = var_names{varname_num};
        temp_var = tab_mean.(curr_varname);
        index_of_negs = temp_var<0;
        temp_var(index_of_negs) = -temp_var(index_of_negs); % all number are positive
        temp_var = arrayfun(@(x) CommaFormat(x),temp_var,'UniformOutput',false);
        temp_var = string(temp_var);
        temp_var(index_of_negs) = cellfun(@(x) ['-' x],temp_var(index_of_negs),'UniformOutput',false); % re-insert negative signs
        tab_mean.(curr_varname) = temp_var;
    end


    %writetable(tab_mean,filename,'Sheet','Region Deaths 2020 2030 mean')


    tab_lower=table(...
        regions_global_list',...
        empty_val_cell, ...
        round(region_excess_deaths_1_2a_mat_lower_2_5), ...
        round(region_excess_deaths_1_2d_mat_lower_2_5), ...
        round(region_excess_deaths_1_3a_mat_lower_2_5), ...
        round(region_excess_deaths_1_3d_mat_lower_2_5), ...
        round(region_excess_deaths_1_4a_mat_lower_2_5), ...
        round(region_excess_deaths_1_4c_mat_lower_2_5) ...
        );


    var_names = {...
        'Region', ...
        'Excess_deaths_1_lower_2_5', ...
        'Excess_deaths_2a_lower_2_5', ...
        'Excess_deaths_2d_lower_2_5', ...
        'Excess_deaths_3a_lower_2_5', ...
        'Excess_deaths_3d_lower_2_5', ...
        'Excess_deaths_4a_lower_2_5', ...
        'Excess_deaths_4c_lower_2_5' ...
        };
    num_varnames = length(var_names);
    tab_lower.Properties.VariableNames = var_names;


    for varname_num = 3:num_varnames
        curr_varname = var_names{varname_num};
        temp_var = tab_lower.(curr_varname);
        index_of_negs = temp_var<0;
        temp_var(index_of_negs) = -temp_var(index_of_negs); % all number are positive
        temp_var = arrayfun(@(x) CommaFormat(x),temp_var,'UniformOutput',false);
        temp_var = string(temp_var);
        temp_var(index_of_negs) = cellfun(@(x) ['-' x],temp_var(index_of_negs),'UniformOutput',false); % re-insert negative signs
        tab_lower.(curr_varname) = temp_var;
    end


    %writetable(tab_lower,filename,'Sheet','Region Deaths 2020 2030 lower 2.5')


    tab_upper=table(...
        regions_global_list',...
        empty_val_cell, ...
        round(region_excess_deaths_1_2a_mat_upper_97_5), ...
        round(region_excess_deaths_1_2d_mat_upper_97_5), ...
        round(region_excess_deaths_1_3a_mat_upper_97_5), ...
        round(region_excess_deaths_1_3d_mat_upper_97_5), ...
        round(region_excess_deaths_1_4a_mat_upper_97_5), ...
        round(region_excess_deaths_1_4c_mat_upper_97_5) ...
        );


    var_names = {...
        'Region', ...
        'Excess_deaths_1_upper_97_5', ...
        'Excess_deaths_2a_upper_97_5', ...
        'Excess_deaths_2d_upper_97_5', ...
        'Excess_deaths_3a_upper_97_5', ...
        'Excess_deaths_3d_upper_97_5', ...
        'Excess_deaths_4a_upper_97_5', ...
        'Excess_deaths_4c_upper_97_5' ...
        };
    num_varnames = length(var_names);
    tab_upper.Properties.VariableNames = var_names;


    for varname_num = 3:num_varnames
        curr_varname = var_names{varname_num};
        temp_var = tab_upper.(curr_varname);
        index_of_negs = temp_var<0;
        temp_var(index_of_negs) = -temp_var(index_of_negs); % all number are positive
        temp_var = arrayfun(@(x) CommaFormat(x),temp_var,'UniformOutput',false);
        temp_var = string(temp_var);
        temp_var(index_of_negs) = cellfun(@(x) ['-' x],temp_var(index_of_negs),'UniformOutput',false); % re-insert negative signs
        tab_upper.(curr_varname) = temp_var;
    end


    %writetable(tab_upper,filename,'Sheet','Region deaths 2020 2030 upper 97.5')


    tab = build_mean_min_max_string(tab_mean,tab_lower,tab_upper,num_regions_global,8);


    var_names = {...
        'Region', ...
        'Excess_deaths_1', ...
        'Excess_deaths_2a', ...
        'Excess_deaths_2d', ...
        'Excess_deaths_3a', ...
        'Excess_deaths_3d', ...
        'Excess_deaths_4a', ...
        'Excess_deaths_4c' ...
        };
    num_varnames = length(var_names);
    tab.Properties.VariableNames = var_names;


    writetable(tab,filename,'Sheet','Region Excess deaths 2020 2030')


end % end of draw_Suppl_Tables_3_4_fun function






function [diff_scenarios_mean_mat,diff_scenarios_lower_mat,diff_scenarios_upper_mat,...
    larger_scenario_mean_mat,larger_scenario_lower_mat,larger_scenario_upper_mat,...
    smaller_scenario_mean_mat,smaller_scenario_lower_mat,smaller_scenario_upper_mat] = ...
    difference_in_deaths_fun(sum_cols,...
    countries_list,cohort_years_vec,num_cohort_years,cohort_years_pos_vec,country_groupings,geographical_unit,geographical_unit_num,results_cell_array,num_stochas_runs,results_field_of_interest,...
    larger_scenario_name,smaller_scenario_name,...
    larger_scenario_num,smaller_scenario_num,...
    diff_scenarios_mean_mat,diff_scenarios_lower_mat,diff_scenarios_upper_mat,...
    larger_scenario_mean_mat,larger_scenario_lower_mat,larger_scenario_upper_mat,...
    smaller_scenario_mean_mat,smaller_scenario_lower_mat,smaller_scenario_upper_mat)

        larger_scenario_results_mat = get_model_output(countries_list,cohort_years_vec,cohort_years_pos_vec,...
            country_groupings,geographical_unit,larger_scenario_name,larger_scenario_num,results_cell_array,num_stochas_runs,results_field_of_interest);
        assert(isequal(size(larger_scenario_results_mat),[num_stochas_runs num_cohort_years]))
        if sum_cols
            larger_scenario_results_mat = sum(larger_scenario_results_mat,2);
        end
        if exist('larger_scenario_mean_mat')==1
            larger_scenario_mean_mat(geographical_unit_num,:) = mean(larger_scenario_results_mat,1);
            larger_scenario_lower_mat(geographical_unit_num,:) = prctile(larger_scenario_results_mat,2.5,1);
            larger_scenario_upper_mat(geographical_unit_num,:) = prctile(larger_scenario_results_mat,97.5,1);
        else
            larger_scenario_mean_mat = NaN(num_stochas_runs,num_cohort_years);
            larger_scenario_lower_mat = NaN(num_stochas_runs,num_cohort_years);
            larger_scenario_upper_mat = NaN(num_stochas_runs,num_cohort_years);
        end

        smaller_scenario_results_mat = get_model_output(countries_list,cohort_years_vec,cohort_years_pos_vec,...
            country_groupings,geographical_unit,smaller_scenario_name,smaller_scenario_num,results_cell_array,num_stochas_runs,results_field_of_interest);
        assert(isequal(size(smaller_scenario_results_mat),[num_stochas_runs num_cohort_years]))
        if sum_cols
            smaller_scenario_results_mat = sum(smaller_scenario_results_mat,2);
        end
        if exist('smaller_scenario_mean_mat')==1
            smaller_scenario_mean_mat(geographical_unit_num,:) = mean(smaller_scenario_results_mat,1);
            smaller_scenario_lower_mat(geographical_unit_num,:) = prctile(smaller_scenario_results_mat,2.5,1);
            smaller_scenario_upper_mat(geographical_unit_num,:) = prctile(smaller_scenario_results_mat,97.5,1);
        else
            smaller_scenario_mean_mat = NaN(num_stochas_runs,num_cohort_years);
            smaller_scenario_lower_mat = NaN(num_stochas_runs,num_cohort_years);
            smaller_scenario_upper_mat = NaN(num_stochas_runs,num_cohort_years);
        end

        diff_mat = larger_scenario_results_mat - smaller_scenario_results_mat;
        if sum_cols
            assert(isequal(size(diff_mat),[num_stochas_runs 1]))
        else
            assert(isequal(size(diff_mat),[num_stochas_runs num_cohort_years]))
        end
        diff_scenarios_mean_mat(geographical_unit_num,:) = mean(diff_mat,1);
        diff_scenarios_lower_mat(geographical_unit_num,:) = prctile(diff_mat,2.5,1);
        diff_scenarios_upper_mat(geographical_unit_num,:) = prctile(diff_mat,97.5,1);

end % end of difference_in_deaths_fun function




function [infant_coverage_2020_2030_mat,BD_coverage_2020_2030_mat,Tot_Pop_2025_mat] = ...
    get_infant_BD_Tot_Pop_2025(num_countries,...
    num_time_steps,...
    country_groupings,geographical_unit,geographical_unit_num,results_cell_array,...
    ListOfScenarios,scenario_name,scenario_num,...
    infant_coverage_2020_2030_mat,BD_coverage_2020_2030_mat,Tot_Pop_2025_mat)

    assert(isequal(size(infant_coverage_2020_2030_mat),[num_countries num_time_steps]))
    assert(isequal(size(BD_coverage_2020_2030_mat),[num_countries num_time_steps]))
    assert(isequal(size(Tot_Pop_2025_mat),[num_countries 1]))

    if strcmp(country_groupings,'regions')
        infant_coverage_2020_2030_mat(geographical_unit_num,:) = NaN(1,2);
        BD_coverage_2020_2030_mat(geographical_unit_num,:) = NaN(1,2);
        Tot_Pop_2025_mat(geographical_unit_num,:) = NaN;
    else
        scenario = ListOfScenarios{scenario_num};
        assert(strcmp(scenario,scenario_name))
        scenario_map = results_cell_array(scenario_name);
        scenario_struct = scenario_map(geographical_unit);

        infant_vacc_rate_vec = scenario_struct.InfantVacc * 100;
        assert(isequal(size(infant_vacc_rate_vec),[1 num_time_steps]))
        infant_coverage_2020_2030_mat(geographical_unit_num,:) = infant_vacc_rate_vec;

        BDvacc_rate_vec = scenario_struct.BirthDoseVacc * 100;
        assert(isequal(size(BDvacc_rate_vec),[1 num_time_steps]))
        BD_coverage_2020_2030_mat(geographical_unit_num,:) = BDvacc_rate_vec;

        Tot_Pop_2025_mat(geographical_unit_num,:) = scenario_struct.Tot_Pop_2025;
    end

    assert(isequal(size(infant_coverage_2020_2030_mat),[num_countries num_time_steps]))
    assert(isequal(size(BD_coverage_2020_2030_mat),[num_countries num_time_steps]))
    assert(isequal(size(Tot_Pop_2025_mat),[num_countries 1]))

end




function [BD_coverage_2020_2030_mat,Tot_Pop_2025_mat,...
    NumSAg_chronic_5_year_olds_mat,Tot_Pop_5_year_olds_mat] = ...
    get_BD_NumSAg_Tot_Pop(num_countries,countries_list,...
    num_year_divisions,calandar_years_vec,num_calandar_years,calandar_years_pos_vec,begin_year_index,end_year_index,...
    country_groupings,geographical_unit,geographical_unit_num,results_cell_array,vaccination_map_default,num_stochas_runs,...
    ListOfScenarios,scenario_name,scenario_num,...
    BD_coverage_2020_2030_mat,Tot_Pop_2025_mat,...
    NumSAg_chronic_5_year_olds_mat,Tot_Pop_5_year_olds_mat)

        assert(isequal(size(BD_coverage_2020_2030_mat),[num_countries 2]))
        assert(isequal(size(NumSAg_chronic_5_year_olds_mat),[num_countries num_stochas_runs num_calandar_years]))
        assert(isequal(size(Tot_Pop_5_year_olds_mat),[num_countries num_stochas_runs num_calandar_years]))

        if strcmp(country_groupings,'regions')
            BD_coverage_2020_2030_mat(geographical_unit_num,:) = NaN(1,2);
            Tot_Pop_2025_mat(geographical_unit_num,:) = NaN;
        else
            scenario = ListOfScenarios{scenario_num};
            assert(strcmp(scenario,scenario_name))
            scenario_map = vaccination_map_default(scenario_name);
            scenario_struct = scenario_map(geographical_unit);
            BDvacc_rate_vec = scenario_struct.BirthDoseVacc(1:num_year_divisions:end) * 100;
            assert(isequal(size(BDvacc_rate_vec),[1 (2100-1979)]))
            BD_coverage_2020_2030_mat(geographical_unit_num,:) = BDvacc_rate_vec([begin_year_index end_year_index]);
            Tot_Pop_2025_mat(geographical_unit_num,:) = scenario_struct.Tot_Pop_2025;
        end

        NumSAg_chronic_5_year_olds_mat(geographical_unit_num,:,:) = get_model_output(countries_list,calandar_years_vec,calandar_years_pos_vec,...
            country_groupings,geographical_unit,scenario_name,scenario_num,results_cell_array,num_stochas_runs,'NumSAg_chronic_1yr_5_year_olds');
        Tot_Pop_5_year_olds_mat(geographical_unit_num,:,:) = get_model_output(countries_list,calandar_years_vec,calandar_years_pos_vec,...
            country_groupings,geographical_unit,scenario_name,scenario_num,results_cell_array,num_stochas_runs,'Tot_Pop_1yr_5_year_olds');

        assert(isequal(size(BD_coverage_2020_2030_mat),[num_countries 2]))
        assert(isequal(size(Tot_Pop_2025_mat),[num_countries 1]))
        assert(isequal(size(NumSAg_chronic_5_year_olds_mat),[num_countries num_stochas_runs num_calandar_years]))
        assert(isequal(size(Tot_Pop_5_year_olds_mat),[num_countries num_stochas_runs num_calandar_years]))

end




function year_of_elim_prev_mat_lower_median_upper = ...
    find_year_of_elim_fun(countries_list,calandar_years_vec,calandar_years_pos_vec,...
    country_groupings,geographical_unit,geographical_unit_num,results_cell_array,num_stochas_runs,threshold_val,...
    scenario_name,scenario_num,...
    year_of_elim_prev_mat_lower_median_upper)

        prev_perc_mat = get_model_output(countries_list,calandar_years_vec,calandar_years_pos_vec,...
            country_groupings,geographical_unit,scenario_name,scenario_num,results_cell_array,num_stochas_runs,'NumSAg_chronic_1yr_5_year_olds','Tot_Pop_1yr_5_year_olds');
        prev_perc_mat = prev_perc_mat * 100;
        year_of_elim_prev_vec = arrayfun(@(xx) ...
            assign_year_of_elim(calandar_years_vec,prev_perc_mat(xx,:),threshold_val), ...
            1:size(prev_perc_mat,1))';
        assert(isequal(size(year_of_elim_prev_vec),[num_stochas_runs 1]))
        year_of_elim_prev_mat_lower_median_upper(geographical_unit_num,2) = median(year_of_elim_prev_vec,1);
        year_of_elim_prev_mat_lower_median_upper(geographical_unit_num,1) = prctile(year_of_elim_prev_vec,2.5,1);
        year_of_elim_prev_mat_lower_median_upper(geographical_unit_num,3) = prctile(year_of_elim_prev_vec,97.5,1);

end




function out_array = get_model_output(ListOfISOs,years_vec,years_pos_vec,...
    country_groupings,grouping_var,scenario,scenario_num,results_object,num_rows,...
    model_result_numerator,model_result_denomenator)

    num_years = length(years_vec);
    assert(num_years==length(years_pos_vec))

    switch country_groupings
        case 'countries'
            assert(ismember(grouping_var,ListOfISOs))
        case 'regions'
            assert(ismember(grouping_var,{'AFRO','EMRO','EURO','PAHO','SEARO','WPRO','Global'}))
    end

    scenario_map = results_object{scenario_num};
    scenario_struct = scenario_map(grouping_var);
    if exist('model_result_denomenator')==1
        scenario_array = scenario_struct.(model_result_numerator) ./ scenario_struct.(model_result_denomenator);
    else
        scenario_array = scenario_struct.(model_result_numerator);
    end
    out_array = scenario_array(:,years_pos_vec);
    assert(isequal(size(out_array),[num_rows num_years]))

end




function annotate_subplot(plot_pos_vec,plot_label,num_plots_across)

    assert(length(plot_pos_vec)==4)
    assert(num_plots_across>0)
    assert(num_plots_across<3)

    if num_plots_across==2
        pos_vec = [plot_pos_vec(1)-0.06 plot_pos_vec(2)+0.05 plot_pos_vec(3:4)];
    else
        assert(num_plots_across==1)
        pos_vec = [plot_pos_vec(1)-0.06*2 plot_pos_vec(2)+0.05*1.5 plot_pos_vec(3:4)];
    end

    a = annotation('textbox','String',plot_label,'Position',pos_vec,'EdgeColor','none');
    a.FontSize = 7;

end




function ans_val = assign_year_of_elim(years_vec,vec_to_be_tested,val_for_testing)

    num_years = length(years_vec);
    assert(length(vec_to_be_tested)==num_years)
    assert(length(val_for_testing)==1)

    year_of_elim_pos = find(vec_to_be_tested<=val_for_testing,1);
    if isempty(year_of_elim_pos)
        ans_val = 2101; % indicates a censored value
    else
        ans_val = years_vec(year_of_elim_pos);
    end

end




function adjust_years_fun(lower_num,median_num,upper_num)

    assert(median_num>=lower_num)
    assert(median_num<=upper_num)

    median_num = num2str(median_num);
    lower_num = num2str(lower_num);
    upper_num = num2str(upper_num);

    if strcmp(median_num,'1980') median_num = 'before 1980'; end
    if strcmp(lower_num,'1980') lower_num = 'before 1980'; end
    if strcmp(upper_num,'1980') upper_num = 'before 1980'; end
    if strcmp(median_num,'2101') median_num = 'after 2100'; end
    if strcmp(lower_num,'2101') lower_num = 'after 2100'; end
    if strcmp(upper_num,'2101') upper_num = 'after 2100'; end

    disp([median_num ' (' lower_num ', ' upper_num ')'])

end




function [commaFormattedString] = CommaFormat(value)
  % Split into integer part and fractional part.
  [integerPart, decimalPart]=strtok(num2str(value),'.'); 
  % Reverse the integer-part string.
  integerPart=integerPart(end:-1:1); 
  % Insert commas every third entry.
  integerPart=[sscanf(integerPart,'%c',[3,inf])' ... 
      repmat(',',ceil(length(integerPart)/3),1)]'; 
  integerPart=integerPart(:)'; 
  % Strip off any trailing commas.
  integerPart=deblank(integerPart(1:(end-1)));
  % Piece the integer part and fractional part back together again.
  commaFormattedString = [integerPart(end:-1:1) decimalPart];
  return; % CommaFormat
end
 

       

function out_tab = build_mean_min_max_string(tab_mean,tab_lower,tab_upper,num_rows,num_cols)

    assert(isequal(size(tab_mean),[num_rows num_cols]))
    assert(isequal(size(tab_lower),[num_rows num_cols]))
    assert(isequal(size(tab_upper),[num_rows num_cols]))
    assert(isstring(tab_mean{1,4}))

    row_index_vec = 1:size(tab_mean,1);
    out_tab = tab_mean;
    for col_num=3:num_cols
        if strcmp(out_tab{2,col_num},'-')
            continue % pass control to next iteration of for or while loop
        else
            out_tab(:,col_num) = arrayfun(@(xxx) ...
                tab_mean{xxx,col_num} + " (" + tab_lower{xxx,col_num} + " to " + tab_upper{xxx,col_num} + ")",row_index_vec,...
                'UniformOutput',false)';
        end
    end
    assert(isequal(size(out_tab),[num_rows num_cols]))

end




function [stochas_regions_cell_array,stochas_countries_cell_array] = make_stochas_regions_countries_cell_arrays(...
    infilepath_stochastic_results,file_string,...
    num_scenarios,scenario_nums_struct)

    [results_regions_scenario_status_quo_stochas,results_countries_scenario_status_quo_stochas] = get_regions_countries_structures(infilepath_stochastic_results,file_string,'1');
    [results_regions_scenario_BD_25_stochas,results_countries_scenario_BD_25_stochas] = get_regions_countries_structures(infilepath_stochastic_results,file_string,'2a');
    [results_regions_scenario_BD_50_stochas,results_countries_scenario_BD_50_stochas] = get_regions_countries_structures(infilepath_stochastic_results,file_string,'2b');
    [results_regions_scenario_BD_75_stochas,results_countries_scenario_BD_75_stochas] = get_regions_countries_structures(infilepath_stochastic_results,file_string,'2c');
    [results_regions_scenario_BD_90_stochas,results_countries_scenario_BD_90_stochas] = get_regions_countries_structures(infilepath_stochastic_results,file_string,'2d');
    [results_regions_scenario_BD_drop_5_stochas,results_countries_scenario_BD_drop_5_stochas] = get_regions_countries_structures(infilepath_stochastic_results,file_string,'3a');
    [results_regions_scenario_BD_drop_10_stochas,results_countries_scenario_BD_drop_10_stochas] = get_regions_countries_structures(infilepath_stochastic_results,file_string,'3b');
    [results_regions_scenario_BD_drop_15_stochas,results_countries_scenario_BD_drop_15_stochas] = get_regions_countries_structures(infilepath_stochastic_results,file_string,'3c');
    [results_regions_scenario_BD_drop_20_stochas,results_countries_scenario_BD_drop_20_stochas] = get_regions_countries_structures(infilepath_stochastic_results,file_string,'3d');
    [results_regions_scenario_2023_2030_stochas,results_countries_scenario_2023_2030_stochas] = get_regions_countries_structures(infilepath_stochastic_results,file_string,'4a');
    [results_regions_scenario_2023_2033_stochas,results_countries_scenario_2023_2033_stochas] = get_regions_countries_structures(infilepath_stochastic_results,file_string,'4b');
    [results_regions_scenario_2025_2040_stochas,results_countries_scenario_2025_2040_stochas] = get_regions_countries_structures(infilepath_stochastic_results,file_string,'4c');


    stochas_regions_cell_array{scenario_nums_struct.status_quo_infant_status_quo_BD_num} = results_regions_scenario_status_quo_stochas;
    stochas_regions_cell_array{scenario_nums_struct.status_quo_infant_BD_to_25_num} = results_regions_scenario_BD_25_stochas;
    stochas_regions_cell_array{scenario_nums_struct.status_quo_infant_BD_to_50_num} = results_regions_scenario_BD_50_stochas;
    stochas_regions_cell_array{scenario_nums_struct.status_quo_infant_BD_to_75_num} = results_regions_scenario_BD_75_stochas;
    stochas_regions_cell_array{scenario_nums_struct.status_quo_infant_BD_to_90_num} = results_regions_scenario_BD_90_stochas;
    stochas_regions_cell_array{scenario_nums_struct.status_quo_infant_BD_drop_5_num} = results_regions_scenario_BD_drop_5_stochas;
    stochas_regions_cell_array{scenario_nums_struct.status_quo_infant_BD_drop_10_num} = results_regions_scenario_BD_drop_10_stochas;
    stochas_regions_cell_array{scenario_nums_struct.status_quo_infant_BD_drop_15_num} = results_regions_scenario_BD_drop_15_stochas;
    stochas_regions_cell_array{scenario_nums_struct.status_quo_infant_BD_drop_20_num} = results_regions_scenario_BD_drop_20_stochas;
    stochas_regions_cell_array{scenario_nums_struct.BD_delayed_2023_2030_num} = results_regions_scenario_2023_2030_stochas;
    stochas_regions_cell_array{scenario_nums_struct.BD_delayed_slowed_2023_2033_num} = results_regions_scenario_2023_2033_stochas;
    stochas_regions_cell_array{scenario_nums_struct.BD_delayed_slowed_2025_2040_num} = results_regions_scenario_2025_2040_stochas;
    assert(length(stochas_regions_cell_array)==num_scenarios)


    stochas_countries_cell_array{scenario_nums_struct.status_quo_infant_status_quo_BD_num} = results_countries_scenario_status_quo_stochas;
    stochas_countries_cell_array{scenario_nums_struct.status_quo_infant_BD_to_25_num} = results_countries_scenario_BD_25_stochas;
    stochas_countries_cell_array{scenario_nums_struct.status_quo_infant_BD_to_50_num} = results_countries_scenario_BD_50_stochas;
    stochas_countries_cell_array{scenario_nums_struct.status_quo_infant_BD_to_75_num} = results_countries_scenario_BD_75_stochas;
    stochas_countries_cell_array{scenario_nums_struct.status_quo_infant_BD_to_90_num} = results_countries_scenario_BD_90_stochas;
    stochas_countries_cell_array{scenario_nums_struct.status_quo_infant_BD_drop_5_num} = results_countries_scenario_BD_drop_5_stochas;
    stochas_countries_cell_array{scenario_nums_struct.status_quo_infant_BD_drop_10_num} = results_countries_scenario_BD_drop_10_stochas;
    stochas_countries_cell_array{scenario_nums_struct.status_quo_infant_BD_drop_15_num} = results_countries_scenario_BD_drop_15_stochas;
    stochas_countries_cell_array{scenario_nums_struct.status_quo_infant_BD_drop_20_num} = results_countries_scenario_BD_drop_20_stochas;
    stochas_countries_cell_array{scenario_nums_struct.BD_delayed_2023_2030_num} = results_countries_scenario_2023_2030_stochas;
    stochas_countries_cell_array{scenario_nums_struct.BD_delayed_slowed_2023_2033_num} = results_countries_scenario_2023_2033_stochas;
    stochas_countries_cell_array{scenario_nums_struct.BD_delayed_slowed_2025_2040_num} = results_countries_scenario_2025_2040_stochas;
    assert(length(stochas_countries_cell_array)==num_scenarios)

end % end function make_stochas_regions_countries_cell_arrays




function [results_regions_scenario,results_countries_scenario] = get_regions_countries_structures(infilepath_stochastic_results,file_string,scenario_num)

    results_regions_countries_scenario = ...
        load([infilepath_stochastic_results 'stochastic_summary_countries_regions_global_' file_string '_scenario_' scenario_num '.mat']);
    % contains countryMap and regionMap
    results_regions_scenario = results_regions_countries_scenario.regionMap;
    results_countries_scenario = results_regions_countries_scenario.countryMap;

end

