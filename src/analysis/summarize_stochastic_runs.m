%%

clear

% this script takes about 01:27:34 hours to summarize 4 sensitivity analyses, each containing 110 countries, 200 runs and 12 scenarios
% this script generates 4 x 12 = 48 results files, mostly about 93 MB each in size
% particles are shuffled within each country before forming regional sums


basedir = fileparts(fileparts(pwd)); % path for folder two levels up from script


sensitivity_analysis_list = {'default','infant_100','treat_medium','treat_high'};
num_sensitivity_analyses = length(sensitivity_analysis_list);


num_runs = 200;


begin_time_run = datetime('now');
disp(['This analysis started at ' datestr(begin_time_run)])


load(fullfile(basedir,'resources','ListOfISOs.mat')) % contains ListOfISOs
load(fullfile(basedir,'resources','WHO_region_map.mat')) % contains WHO_region_map
load(fullfile(basedir,'resources','particle_index_countries.mat')) % contains particle_index_countries_mat


ListOfAllScenarios = {...
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
ListOfAllScenarioNums = {...
    '1',...
    '2a','2b','2c','2d',...
    '3a','3b','3c','3d',...
    '4a','4b','4c'};
num_all_scenarios = length(ListOfAllScenarios);
scenarios_map = containers.Map(ListOfAllScenarioNums, ListOfAllScenarios);
assert(length(keys(scenarios_map))==12)


num_countries = length(ListOfISOs);
assert(num_countries==110)
disp([num2str(num_countries),' countries'])


regions_list = {'AFRO','EMRO','EURO','PAHO','SEARO','WPRO'};
regions_list{end+1} = 'Global';
num_regions = length(regions_list);
assert(num_regions==6+1) % 6 WHO regions and global


num_years = 2100 - 1980 + 1;
num_year_divisions = 10;


% only keeping numbers from 1 to num_runs from particle_index_countries_mat
assert(isequal(size(particle_index_countries_mat),[200 num_countries]))
nums_index = particle_index_countries_mat<=num_runs;
particle_index_countries_mod_mat = reshape(particle_index_countries_mat(nums_index),num_runs,num_countries);
assert(isequal(size(particle_index_countries_mod_mat),[num_runs num_countries]))
if num_runs==200
    assert(isequal(particle_index_countries_mod_mat,particle_index_countries_mat))
end


sensitivity_analysis_hours_vec = repmat(duration(0,0,0),1,num_sensitivity_analyses);


for sensitivity_analysis_num=1:num_sensitivity_analyses

    sensitivity_analysis = sensitivity_analysis_list{sensitivity_analysis_num};
    disp(['Sensitivity analysis: ' sensitivity_analysis])

    begin_time_sensitivity_analysis = datetime('now');
    disp(['The "' sensitivity_analysis '" sensitivity analysis started at ' datestr(begin_time_sensitivity_analysis)])

    for scenario_num = 1:num_all_scenarios

        scenario = ListOfAllScenarios{scenario_num};
        scenario_num_string = ListOfAllScenarioNums{scenario_num};
        disp(['Scenario: ' scenario])

        begin_time_scenario = datetime('now');
        disp(['The "' scenario '" scenario started at ' datestr(begin_time_scenario)])

        summarize_stochastic_runs_sens_scenario(sensitivity_analysis,...
            num_all_scenarios,scenario,scenario_num_string,...
            num_runs,...
            WHO_region_map,regions_list,num_regions,...
            ListOfISOs,num_countries,...
            num_years,...
            particle_index_countries_mod_mat,...
            basedir)

        num_scenarios_left = num_all_scenarios - scenario_num;
        if num_scenarios_left>0        
            disp(['There are ' num2str(num_scenarios_left) ' scenarios left in this sensitivity analysis (' sensitivity_analysis ').'])
        end

    end

    time_taken_for_sensitivity_analysis = datetime('now') - begin_time_sensitivity_analysis;
    disp(['The duration of this sensitivity analysis (' sensitivity_analysis ') was ' char(time_taken_for_sensitivity_analysis) ' hours.\n\n'])
    sensitivity_analysis_hours_vec(sensitivity_analysis_num) = time_taken_for_sensitivity_analysis;
    assert(all(sensitivity_analysis_hours_vec(1:sensitivity_analysis_num)>0))
    average_time_per_sensitivity_analysis = mean(sensitivity_analysis_hours_vec(1:sensitivity_analysis_num));
    min_time_per_sensitivity_analysis = min(sensitivity_analysis_hours_vec(1:sensitivity_analysis_num));
    max_time_per_sensitivity_analysis = max(sensitivity_analysis_hours_vec(1:sensitivity_analysis_num));
    num_sensitivity_analyses_left = num_sensitivity_analyses - sensitivity_analysis_num;
    mean_time_left = num_sensitivity_analyses_left * average_time_per_sensitivity_analysis;
    min_time_left = num_sensitivity_analyses_left * min_time_per_sensitivity_analysis;
    max_time_left = num_sensitivity_analyses_left * max_time_per_sensitivity_analysis;
    if num_sensitivity_analyses_left>0
        disp(['There are ' num2str(num_sensitivity_analyses_left) ' sensitivity analyses left, which will take about ' char(mean_time_left) ' hours.'])
    end

end


end_time_run = datetime('now');
disp(end_time_run)
time_taken_for_run = end_time_run - begin_time_run;
disp(['The duration of this analysis was ' char(time_taken_for_run) ' hours.\n\n'])




function summarize_stochastic_runs_sens_scenario(sensitivity_analysis,...
    num_all_scenarios,scenario,scenario_num_string,...
    num_runs,...
    WHO_region_map,regions_list,num_regions,...
    ListOfISOs,num_countries,...
    num_years,...
    particle_index_countries_mat,...
    basedir)


    num_runs_string = num2str(num_runs);


    begin_string = 'results_countries_';
    end_string = '.mat';
    particles_str = 'stochastic_run_2';    
    infilename = [begin_string sensitivity_analysis '_' particles_str end_string]; % results_countries_stochastic_run_1.mat
    load(fullfile(basedir,'outputs',infilename)); % contains outMap
    assert(length(outMap)==num_all_scenarios) % 12 scenarios
    Runs_map_run_2_scenario_ISO = outMap(scenario); % sensitivity analysis, treatment, scenario, run_num results
    Runs_map_run_2_scenario_ISO = Runs_map_run_2_scenario_ISO(ListOfISOs{1}); % sensitivity analysis, treatment, scenario, run_num AFG's results
    fields_list = sort(fields(Runs_map_run_2_scenario_ISO))';
    num_fields = length(fields_list);
    assert(num_fields==21-4)
    fields_list_without_Time = setdiff(fields_list,{'Time'});
    num_fields_without_Time = length(fields_list_without_Time);
    assert(num_fields-1==num_fields_without_Time)
    fields_out_larger = {'Incid_chronic_all_1yr_approx',...
        'NumSAg_chronic_1yr_5_year_olds','Tot_Pop_1yr_5_year_olds',...
        'Incid_Deaths_1yr_approx','Incid_Deaths_1yr_approx_birth_cohorts'}; % keep all 200 stochastic runs rather that only the summaries
    num_fields_out_larger = length(fields_out_larger);
    assert(all(ismember(fields_out_larger,fields_list)))
    fields_out_larger_1_2d = {'Prev_TDF_treat_1yr','Prev_treatment_eligible_1yr','DALYPerYear'}; % keep all 200 stochastic runs rather that only the summaries
    num_fields_out_larger_1_2d = length(fields_out_larger_1_2d);
    assert(all(ismember(fields_out_larger_1_2d,fields_list)))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.Time),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.Tot_Pop_1yr),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.Incid_chronic_all_1yr_approx),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.NumSAg_1yr),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.NumSAg_chronic_1yr),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.Prev_treatment_eligible_1yr),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.Prev_TDF_treat_1yr),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.Incid_Deaths_1yr_approx),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.DALYPerYear),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.Tot_Pop_1yr_5_year_olds),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.Tot_Pop_1yr_under_5_year_olds),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.NumSAg_1yr_5_year_olds),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.NumSAg_1yr_under_5_year_olds),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.NumSAg_chronic_1yr_5_year_olds),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.NumSAg_chronic_1yr_under_5_year_olds),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.num_births_1yr),[1 num_years]))
    assert(isequal(size(Runs_map_run_2_scenario_ISO.Incid_Deaths_1yr_approx_birth_cohorts),[1 num_years]))
    for nn=1:num_fields
        curr_field = fields_list{nn};
        assert(all(all(all(Runs_map_run_2_scenario_ISO.(curr_field)>=0))))
    end

    
    countryMap = containers.Map; 
    % will contain the countries for this scenario


    %disp('Making country matrices.')


    for run_num=1:num_runs

        stochas_run_str = num2str(run_num);
        disp(run_num)

        particles_str = ['stochastic_run_', stochas_run_str];    
        infilename = [begin_string sensitivity_analysis '_' particles_str end_string]; % results_countries_stochastic_run_1.mat
        load(fullfile(basedir,'outputs',infilename)); % contains outMap
        assert(length(outMap)==num_all_scenarios) % 12 scenarios
        countryMap_source = outMap(scenario);
        assert(length(countryMap_source)==num_countries)

        for country_num=1:num_countries

            ISO = ListOfISOs{country_num};
        
            country_results = countryMap_source(ISO);
            if run_num==1
                country_results = rmfield(country_results,'country_name');
                assert(strcmp(country_results.scenario,scenario))
                country_results = rmfield(country_results,'scenario');
                country_results = rmfield(country_results,'InfantVacc');
                country_results = rmfield(country_results,'BirthDoseVacc');
            end
            assert(isequal(sort(fieldnames(country_results))',fields_list))

            country_results_for_summary_statistics = country_results;
            country_results_for_summary_statistics = rmfield(country_results_for_summary_statistics,'Time');
            assert(isequal(sort(fields(country_results_for_summary_statistics))',fields_list_without_Time))

            country.Time = country_results.Time;
        
            for fieldname_num=1:num_fields_without_Time
                % for each fieldname, construct matrix of model output vectors

                curr_fieldname = fields_list_without_Time{fieldname_num};
                if run_num==1
                    tmp_mat = zeros(num_runs,num_years);
                else
                    country = countryMap(ISO);
                    tmp_mat = country.(curr_fieldname);
                end
                tmp_mat(run_num,:) = country_results_for_summary_statistics.(curr_fieldname);
                country.(curr_fieldname) = tmp_mat;
        
                countryMap(ISO) = country;
                % country is a structure containing a Time vector, as well a matrix for each model output 

            end % end fieldname_num for loop

            assert(length(fieldnames(country))==num_fields)
            assert(isequal(size(country.Time),[1 num_years]))
            assert(isequal(size(country.Incid_chronic_all_1yr_approx),[num_runs num_years]))

            clear 'country' 'country_results' 'country_results_for_summary_statistics'...
                'curr_fieldname' 'fieldname_num' 'tmp_mat'

        end % end country_num for loop
    
        assert(length(countryMap)==num_countries)

    end % end run_num for loop


    %disp('Calculating summary statistics from country matrices.')


    for country_num=1:num_countries

        ISO = ListOfISOs{country_num};

        assert(isequal(size(particle_index_countries_mat),[num_runs num_countries]))
        particle_indices_country_vec = particle_index_countries_mat(:,country_num);
        assert(isequal(size(particle_indices_country_vec),[num_runs 1]))
    
        country_results = countryMap(ISO);
        assert(isequal(sort(fieldnames(country_results))',fields_list))

        country_results_for_summary_statistics = country_results;
        country_results_for_summary_statistics = rmfield(country_results_for_summary_statistics,'Time');
        assert(isequal(sort(fields(country_results_for_summary_statistics))',fields_list_without_Time))

        country.Time = country_results.Time;
    
        for fieldname_num=1:num_fields_without_Time
            curr_fieldname = fields_list_without_Time{fieldname_num};

            if ismember(curr_fieldname,fields_out_larger)
                tmp_mat_country = country_results_for_summary_statistics.(curr_fieldname);
                country.(curr_fieldname) = tmp_mat_country(particle_indices_country_vec,:); % re-order paticles for country before forming regional total
                assert(isequal(size(country.(curr_fieldname)),[num_runs num_years]))
            end
            if ismember(scenario_num_string,{'1','2d'}) && ismember(curr_fieldname,fields_out_larger_1_2d)
                tmp_mat_country = country_results_for_summary_statistics.(curr_fieldname);
                country.(curr_fieldname) = tmp_mat_country(particle_indices_country_vec,:); % re-order paticles for country before forming regional total
                assert(isequal(size(country.(curr_fieldname)),[num_runs num_years]))
            end
    
            countryMap(ISO) = country;

        end % end fieldname_num for loop

        if ismember(scenario_num_string,{'2a','2b','2c','3a','3b','3c','3d','4a','4b','4c'})
            assert(length(fieldnames(country))==1+num_fields_out_larger) % Time and fields_out_larger
        elseif ismember(scenario_num_string,{'1','2d'})
            assert(length(fieldnames(country))==1+num_fields_out_larger+num_fields_out_larger_1_2d) % Time, fields_out_larger and fields_out_larger_1_2d
        else
            error('Unexpected scenario!')
        end

        clear 'country' 'country_results' 'country_results_for_summary_statistics'...
            'curr_fieldname' 'fieldname_num'

    end % end country_num for loop
    
    assert(length(countryMap)==num_countries)


    %disp('Calculating summary statistics for regions.')


    regionMap = containers.Map; 
    % will contain the regions for this scenario

    for country_num=1:num_countries

        ISO = ListOfISOs{country_num};
        region_WHO = WHO_region_map(ISO);
    
        country_results_in = countryMap(ISO);
        country_results_in_fieldnames = sort(fieldnames(country_results_in));
        num_fieldnames_countryMap = length(country_results_in_fieldnames);
        if ismember(scenario_num_string,{'2a','2b','2c','3a','3b','3c','3d','4a','4b','4c'})
            assert(num_fieldnames_countryMap==1+num_fields_out_larger) % Time and fields_out_larger
            concatenated_fieldnames = sort(union({'Time'}, fields_out_larger))';
            assert(isequal(size(concatenated_fieldnames),[(1+num_fields_out_larger) 1]))
            assert(isequal(country_results_in_fieldnames,concatenated_fieldnames))
        elseif ismember(scenario_num_string,{'1','2d'})
            assert(num_fieldnames_countryMap==1+num_fields_out_larger+num_fields_out_larger_1_2d) 
            % Time, fields_out_larger and fields_out_larger_1_2d
            concatenated_fieldnames = sort(union({'Time'}, union(fields_out_larger, fields_out_larger_1_2d)))';
            assert(isequal(size(concatenated_fieldnames),[(1+num_fields_out_larger+num_fields_out_larger_1_2d) 1]))
            assert(isequal(country_results_in_fieldnames,concatenated_fieldnames))
        else
            error('Unexpected scenario!')
        end
    

        if ~ismember(region_WHO,keys(regionMap))
            country_results_in.num_countries = 1;
            regionMap(region_WHO) = country_results_in;
            clear country_results_in
            continue % goes on to next country
        end

       
        country_results_for_summary_statistics = country_results_in;
        country_results_for_summary_statistics = rmfield(country_results_for_summary_statistics,'Time');

        region = regionMap(region_WHO);
        assert(length(fieldnames(region))==num_fieldnames_countryMap+1) % num_countries added
        for fieldname_num_extra=1:num_fields_out_larger
            curr_fieldname_extra = fields_out_larger{fieldname_num_extra};
            region.(curr_fieldname_extra) = country_results_for_summary_statistics.(curr_fieldname_extra) + region.(curr_fieldname_extra);
        end
        if ismember(scenario_num_string,{'1','2d'})
            region.Prev_TDF_treat_1yr = country_results_for_summary_statistics.Prev_TDF_treat_1yr + region.Prev_TDF_treat_1yr;
            region.Prev_treatment_eligible_1yr = country_results_for_summary_statistics.Prev_treatment_eligible_1yr + region.Prev_treatment_eligible_1yr;
            region.DALYPerYear = country_results_for_summary_statistics.DALYPerYear + region.DALYPerYear;
        end
        region.num_countries = region.num_countries + 1;
        if ismember(scenario_num_string,{'2a','2b','2c','3a','3b','3c','3d','4a','4b','4c'})
            assert(length(fieldnames(region))==1+1+num_fields_out_larger) % Time, num_countries and fields_out_larger
        elseif ismember(scenario_num_string,{'1','2d'})
            assert(length(fieldnames(region))==1+1+num_fields_out_larger+num_fields_out_larger_1_2d) 
            % Time, num_countries, fields_out_larger and fields_out_larger_1_2d
        else
            error('Unexpected scenario!')
        end
        regionMap(region_WHO) = region;

        clear 'region' 'country_results_in' 'country_results_for_summary_statistics'...
            'curr_fieldname' 'fieldname_num' 'num_fieldnames_countryMap'

    end % end country_num for loop

    assert(regionMap('AFRO').num_countries==42)
    assert(regionMap('EMRO').num_countries==13)
    assert(regionMap('EURO').num_countries==14)
    assert(regionMap('PAHO').num_countries==15)
    assert(regionMap('SEARO').num_countries==10)
    assert(regionMap('WPRO').num_countries==16)


    %disp('Calculating global summary statistics.')
 

    for region_num=1:(num_regions-1) % leaving out 'Global'

        region_WHO = regions_list{region_num};
    
        region_model_results_in = regionMap(region_WHO);
        region_model_results_in_fieldnames = sort(fieldnames(region_model_results_in));
        num_fieldnames_regionMap = length(region_model_results_in_fieldnames);
        if ismember(scenario_num_string,{'2a','2b','2c','3a','3b','3c','3d','4a','4b','4c'})
            assert(num_fieldnames_regionMap==1+1+num_fields_out_larger) % Time, num_countries and fields_out_larger
            concatenated_fieldnames = sort(union({'Time','num_countries'}, fields_out_larger))';
            assert(isequal(size(concatenated_fieldnames),[(1+1+num_fields_out_larger) 1]))
            assert(isequal(region_model_results_in_fieldnames,concatenated_fieldnames))
        elseif ismember(scenario_num_string,{'1','2d'})
            assert(num_fieldnames_regionMap==1+1+num_fields_out_larger+num_fields_out_larger_1_2d) 
            % Time, num_countries, fields_out_larger and fields_out_larger_1_2d
            concatenated_fieldnames = sort(union({'Time','num_countries'}, union(fields_out_larger, fields_out_larger_1_2d)))';
            assert(isequal(size(concatenated_fieldnames),[(1+1+num_fields_out_larger+num_fields_out_larger_1_2d) 1]))
            assert(isequal(region_model_results_in_fieldnames,concatenated_fieldnames))
        else
            error('Unexpected scenario!')
        end
    

        if ~ismember('Global',keys(regionMap))
            regionMap('Global') = region_model_results_in;
            clear region_model_results_in
            continue % goes on to next country
        end

       
        region_results_for_summary_statistics = region_model_results_in;
        region_results_for_summary_statistics = rmfield(region_results_for_summary_statistics,'Time');

        global_struct = regionMap('Global');
        assert(length(fieldnames(global_struct))==num_fieldnames_regionMap) % Time and num_countries
        global_struct.num_countries = region_results_for_summary_statistics.num_countries + global_struct.num_countries;
        for fieldname_num_extra=1:num_fields_out_larger
            curr_fieldname_extra = fields_out_larger{fieldname_num_extra};
            global_struct.(curr_fieldname_extra) = region_results_for_summary_statistics.(curr_fieldname_extra) + global_struct.(curr_fieldname_extra);
        end
        if ismember(scenario_num_string,{'1','2d'})
            global_struct.Prev_TDF_treat_1yr = region_results_for_summary_statistics.Prev_TDF_treat_1yr + global_struct.Prev_TDF_treat_1yr;
            global_struct.Prev_treatment_eligible_1yr = region_results_for_summary_statistics.Prev_treatment_eligible_1yr + global_struct.Prev_treatment_eligible_1yr;
            global_struct.DALYPerYear = region_results_for_summary_statistics.DALYPerYear + global_struct.DALYPerYear;
        end
        if ismember(scenario_num_string,{'2a','2b','2c','3a','3b','3c','3d','4a','4b','4c'})
            assert(length(fieldnames(global_struct))==1+1+num_fields_out_larger) % Time, num_countries and fields_out_larger
        elseif ismember(scenario_num_string,{'1','2d'})
            assert(length(fieldnames(global_struct))==1+1+num_fields_out_larger+num_fields_out_larger_1_2d) 
            % Time, num_countries, fields_out_larger and fields_out_larger_1_2d
        else
            error('Unexpected scenario!')
        end
        regionMap('Global') = global_struct;

        clear 'global_struct' 'region_model_results_in' 'region_results_for_summary_statistics'...
            'curr_fieldname' 'fieldname_num' 'num_fieldnames_regionMap'

    end % end region_num for loop
       
    assert(length(regionMap)==num_regions)
    assert(regionMap('Global').num_countries==num_countries)

    save(...
        fullfile(basedir,'outputs',['stochastic_summary_countries_regions_global_' sensitivity_analysis '_scenario_', scenario_num_string '.mat']),...
        'countryMap','regionMap')


end % end summarize_stochastic_runs_sens_scenario function

