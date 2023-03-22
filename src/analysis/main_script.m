function main_script(do_geometric_median,save_results,replacement_ListOfISOs,BD_scaleup_cov,PAP_scaleup_end_year,...
    parameter_to_vary,limit_to_AFRO_BD_0,replacement_list_of_strategies)

% To generate results for Fig. 1A, 2A and 2B and Table 3 and Supplementary Table 5
%     main_script(false,true)
% To generate results for Fig. 1B,
%     main_script(true,true)
%     For a = 1, this script takes 01:18:59 hours on desktop computer for 110 countries, 1 particle per country and 6 strategies
% To generate results for Fig. 3,
%     main_script(true,true,[],1,2025)
% To generate results for Fig. 4,
%     main_script(true,true,[],[],[],'pr_VerticalTransmission_HbSAgHighVL_PAP',true,{'1','3','4'})



    begin_time = datetime('now');
    disp(begin_time)

    tic



    aaa = 'a1';
    disp(['Assumption set: ' aaa])


    disp(char(join(['do_geometric_median:' string(do_geometric_median)])))


    if do_geometric_median
        num_stochas_runs_to_run = 1;
        suppress_screen_output = false;
    else
        num_stochas_runs_to_run = 200;
        suppress_screen_output = true;
    end
    disp(['number of runs: ' num2str(num_stochas_runs_to_run)])


    treatment_string = 'no_treat';
    do_treatment = false;
    disp(['treatment: ' treatment_string])



    disp('Importing data structures')


    currentFolder = pwd;
    pat = fullfile('icl-hbv','src','analysis');
    if ~endsWith(currentFolder,pat)
        warning(['Please run this script from within the folder ' pat])
    end
    basedir = fileparts(fileparts(currentFolder)); % path for folder two levels up from script
    addpath(fullfile(basedir,'src','analysis')) % path for function country_level_analyses
    addpath(fullfile(basedir,'src','model')) % path for function HBVmodel

    folder_path_out = fullfile(basedir,'outputs');

    load(fullfile(basedir,'resources','WHO_region_map.mat')) % contains WHO_region_map
    load(fullfile(basedir,'resources','ListOfISOs.mat')) % contains ListOfISOs
    load(fullfile(basedir,'resources','vaccination_coverages.mat')) % contains BD_table and HepB3_table
    load(fullfile(basedir,'resources','reference_prevalences.mat')) % contains country_s_e_HCCdeaths_map
    load(fullfile(basedir,'resources','params_map.mat')) % contains params_map and dwvec
    load(fullfile(basedir,'resources','stochastic_parameters_mat.mat')) % contains stochas_params_mat
    load(fullfile(basedir,'resources','treatment_2016_map.mat')) % contains num_in_treatment_2016_map and pop_size_HBsAg_treatment_map
    treatment_rates_map = load(fullfile(basedir,'resources','treatment_rates_map.mat')); % contains treatment_rates_map_geo_med and treatment_rates_map_stochas
    if do_geometric_median % if do_geometric_median is true
        treatment_rates_map = treatment_rates_map.treatment_rates_map_geo_med;
    else
        treatment_rates_map = treatment_rates_map.treatment_rates_map_stochas;
    end

    disp('Finished importing data structures')



    if ~save_results
        disp('Checking particle values.')
    end


    if exist('replacement_ListOfISOs')==1 % i.e. is a variable in the workspace
        if ~isempty(replacement_ListOfISOs)
            disp('List of countries:')
            disp(replacement_ListOfISOs)
        end
    else
        replacement_ListOfISOs = []; % make an empty object
    end


    ListOfAllISOs = ListOfISOs;
    num_all_countries = length(ListOfAllISOs);
    assert(num_all_countries==110)
    if ~isempty(replacement_ListOfISOs)
        ListOfISOs = replacement_ListOfISOs;
        % do all ISOs specified in replacement_ListOfISOs that are listed under this scenario
    end
    assert(all(ismember(ListOfISOs,ListOfAllISOs)))
    num_countries = length(ListOfISOs);
    assert(num_countries<=num_all_countries)
    if ~suppress_screen_output
        if num_countries>1
            disp([num2str(num_countries),' countries'])
        else
            assert(num_countries>0)
            disp([num2str(num_countries),' country'])
        end
    end


    % variables for importing values for the 7 fitted parameters for a particular country-particle combination
    country_start_cols = 2:8:(1+num_all_countries*8); % there are 1 + num_all_countries*8 + 3 columns in stochas_params_mat


    if exist('BD_scaleup_cov')~=1 % i.e. is not a variable in the workspace
        BD_scaleup_cov = 0.9;
    else
        if isempty(BD_scaleup_cov)
            BD_scaleup_cov = 0.9;
        end
    end
    assert(BD_scaleup_cov>=0)
    assert(BD_scaleup_cov<=1)

    TScaleup_BirthDoseVaccineIntv_start = 2022;
    TScaleup_BirthDoseVaccineIntv_finish = 2024;
    BD_scaleup_cov_str = num2str(BD_scaleup_cov * 100);


    if exist('PAP_scaleup_end_year')~=1 % i.e. is not a variable in the workspace
        PAP_scaleup_end_year = 2024;
    else
        if isempty(PAP_scaleup_end_year)
            PAP_scaleup_end_year = 2024;
        end
    end
    assert(PAP_scaleup_end_year>2019)
    assert(PAP_scaleup_end_year<=2100)

    if ~suppress_screen_output
        disp(['Birth dose coverage scaled up to ' BD_scaleup_cov_str '% between ' num2str(TScaleup_BirthDoseVaccineIntv_start) ' and ' num2str(TScaleup_BirthDoseVaccineIntv_finish) '.'])
        disp(['PAP coverage scaled up to 100% between ' num2str(PAP_scaleup_end_year-1) ' and ' num2str(PAP_scaleup_end_year) '.'])
    end


    if exist('parameter_to_vary')~=1 % i.e. is not a variable in the workspace
        parameter_to_vary = []; % make an empty object
    end


    if exist('limit_to_AFRO_BD_0')~=1 % i.e. is not a variable in the workspace
        limit_to_AFRO_BD_0 = false; % analyse all countries rather than only AFRO countries with no birth dose coverage
    end


    if exist('replacement_list_of_strategies')==1 % i.e. is a variable in the workspace
        if ~isempty(replacement_list_of_strategies)
            assert(~isempty(parameter_to_vary))
            assert(ismember('1',replacement_list_of_strategies)) % ensure that Strategy 1 is one of the strategies being done
        end
    else
        replacement_list_of_strategies = []; % make an empty object
    end


    if save_results
        num_strategies = 6;
    else
        num_strategies = 1;
    end
    disp(['number of strategies in total: ' num2str(num_strategies)])
    if ~isempty(replacement_list_of_strategies)
        disp('Strategies to analyse:')
        disp(replacement_list_of_strategies)
    end


    start_year_simulations = 1890;
    start_year_vacc_coverage = 1980;
    if save_results
        end_year_simulations = 2101; % must go up to year 2101 to get incidence readings up to year 2100
        first_year_results = 1980;
        last_year_results = end_year_simulations - 1; % 2100
        num_years_results = last_year_results - first_year_results + 1; % 121
        T0 = end_year_simulations - start_year_simulations; % 211
    else
        end_year_simulations = 1900;
        first_year_results = start_year_simulations; % 1890
        last_year_results = end_year_simulations - 1; % 1899
        num_years_results = last_year_results - first_year_results + 1; % 1899 - 1890 + 1 = 10; T0: 10
        T0 = end_year_simulations - start_year_simulations; % 10; need to run the model to ensure that everything is working; do not need any results out of the model
    end
    num_year_divisions = 10;
    assert(rem(log10(num_year_divisions),1)==0) % ensure that num_year_divisions is a multiple of 10
    % log10(1) = 0, log10(10) = 1, log10(100) = 2, etc.
    % rem(x,1) gives 0 if x is an integer
    dt = 1/num_year_divisions; % time step

    ages = 0:dt:(100 - dt); % The age group at the beginning of each compartment; 1 x 1000 double; [0 0.1 0.2 ... 99.8 99.9]
    num_age_steps = length(ages); % Number of age-categories (assuming max age is 100 years (i.e. oldest person we see at the start of an iteration is aged 100-dt)
    i5y = 6; % index for 5-year-olds

    theta = 0.01; % proportion of persons newly infected that develop severe symptoms during acute infection
    CFR_Acute = 0.05; % proportion of people with severe acute infection that die from severe acute infection
    rate_6months = 2; % the rate necessary to move everyone out of States 14 and 15 into State 2 in 6 months

    ECofactor = 15; % Multiple for rate of transmission for horizontal transmission with HbEAg vs HbSAg
    % epsilon_2 and epsilon_3 in Table S2, p. 6 of the appendix to the publication


    Snames = {...
        'Susceptible', ...                 % 1
        'HBV: Immune Tolerant', ...        % 2
        'HBV: Immune Reactive', ...        % 3
        'HBV: Asymptomatic Carrier', ...   % 4
        'HBV: Chronic Hep B', ...          % 5
        'HBV: Comp Cirrhosis', ...         % 6
        'HBV: Decomp Cirrhosis', ...       % 7
        'HBV: Liver Cancer', ...           % 8
        'HBV: Immune (Rec. or vacc.)', ... % 9
        'HBV: TDF-Treatment', ...          % 10
        'Prematurely dead due to HBV', ... % 11
        '3TC-Treatment', ...               % 12
        'Failed 3TC-Treatment', ...        % 13
        'Non-severe acute', ...            % 14
        'Severe acute' ...                 % 15
        }; % 1 x 15 cell array
    num_states = length(Snames);
    assert(num_states==15) % 15 disease states


    % risk of chronic carriage due to horizontal and vertical transmission
    risk_ages_edmunds_func = @(age) exp(-0.645*age^(0.455));
    p_ChronicCarriage = 0.885 * ones(1, num_age_steps, 2, 2);
    % does not differ over disease stages
    % 4 dimensions necessary so that it can be multiplied by X
    risk_p_ChronicCarriage = arrayfun(@(xx) risk_ages_edmunds_func(xx), ages(:,ages>0.5));
    risk_p_ChronicCarriage = repmat(risk_p_ChronicCarriage,1,1,2,2);
    p_ChronicCarriage(:, ages>0.5, :, :) = risk_p_ChronicCarriage;


    stochas_run_hours_vec = -99 * ones(1,num_stochas_runs_to_run);


    for stochas_run_num=1:num_stochas_runs_to_run


        stochas_run_str = num2str(stochas_run_num);
        disp(stochas_run_str)

    
        sub_begin_time = datetime('now');
        %disp(sub_begin_time)


        country_level_analyses(...
            suppress_screen_output,aaa,do_geometric_median,...
            stochas_run_num,stochas_run_str,...
            do_treatment,treatment_string,...
            num_strategies,save_results,folder_path_out,...
            start_year_simulations,end_year_simulations,T0,start_year_vacc_coverage,first_year_results,last_year_results,num_years_results,num_year_divisions,dt,...
            ages,num_age_steps,i5y,...
            theta,CFR_Acute,rate_6months,ECofactor,p_ChronicCarriage,...
            num_states,...
            ListOfISOs,num_countries,country_start_cols,...
            BD_scaleup_cov,BD_scaleup_cov_str,TScaleup_BirthDoseVaccineIntv_start,TScaleup_BirthDoseVaccineIntv_finish,...
            PAP_scaleup_end_year,parameter_to_vary,limit_to_AFRO_BD_0,replacement_list_of_strategies,...
            WHO_region_map,ListOfAllISOs,BD_table,HepB3_table,...
            country_s_e_HCCdeaths_map,dwvec,params_map,stochas_params_mat,...
            num_in_treatment_2016_map,pop_size_HBsAg_treatment_map,treatment_rates_map)


        if ~do_geometric_median
            time_taken_for_stochas_run = hours(datetime('now') - sub_begin_time);
            stochas_run_hours_vec(stochas_run_num) = time_taken_for_stochas_run;
            assert(all(stochas_run_hours_vec(1:stochas_run_num)>0))
            average_time_per_stochas_run = mean(stochas_run_hours_vec(1:stochas_run_num));
            min_time_per_stochas_run = min(stochas_run_hours_vec(1:stochas_run_num));
            max_time_per_stochas_run = max(stochas_run_hours_vec(1:stochas_run_num));
            num_stochas_runs_left = num_stochas_runs_to_run - stochas_run_num;
            mean_time_left = num_stochas_runs_left * average_time_per_stochas_run;
            min_time_left = num_stochas_runs_left * min_time_per_stochas_run;
            max_time_left = num_stochas_runs_left * max_time_per_stochas_run;
            if num_stochas_runs_left>10
                if rem(num_stochas_runs_left, 10) == 0 % only show for every 10th country
                    disp([...
                        'There are ' ...
                        num2str(num_stochas_runs_left) ...
                        ' particles left, which will take about ' num2str(mean_time_left) ' (' num2str(min_time_left) ' to ' num2str(max_time_left) ') hours.'...
                        ])
                end
            else
                disp([...
                    'There are ' ...
                    num2str(num_stochas_runs_left) ...
                    ' particles left, which will take about ' num2str(mean_time_left) ' (' num2str(min_time_left) ' to ' num2str(max_time_left) ') hours.'...
                    ])
            end
        end


    end % end for stochas_run_num loop



    toc



    disp('Run finished!')

    end_time = datetime('now');
    disp(end_time)
    disp(['The duration of this run was ' char(end_time - begin_time) '\n\n'])


end % end function main_script


