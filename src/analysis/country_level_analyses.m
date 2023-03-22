function country_level_analyses(...
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

% For a = 1, this script takes 01:16:25 hours on desktop computer for 110 countries, 1 particle per country and 6 strategies


    if nargin < 1
       error('No input')
    end


    if ~suppress_screen_output
        begin_time = datetime('now');
        disp(begin_time)

        tic
    end


    if ~isempty(parameter_to_vary)
        assert(do_geometric_median)
        assert(save_results)
        assert(strcmp(parameter_to_vary,'pr_VerticalTransmission_HbSAgHighVL_PAP'))
        parameter_to_vary_values_vec = [0.05 0.06 0.1:0.1:1];
        parameter_to_vary_num_vals = length(parameter_to_vary_values_vec);
        other_parameter_to_vary = strrep(parameter_to_vary,'HbSAg','HbEAg');
        assert(~isempty(other_parameter_to_vary))

        parameter_to_vary_short = replace(parameter_to_vary,'pr_VerticalTransmission_','');
        parameter_to_vary_short = replace(parameter_to_vary_short,'HbSAg','');
        parameter_to_vary_short = replace(parameter_to_vary_short,'HbEAg','');
        parameter_to_vary_short = ['_' parameter_to_vary_short];
    else
        parameter_to_vary_num_vals = 1; % each model parameter has a single value

        parameter_to_vary_short = '';
    end


    % Select one set of assumptions
    if ~suppress_screen_output
        if ~isempty(parameter_to_vary)
            disp(['Model parameters to vary: ' parameter_to_vary ' and ' other_parameter_to_vary])
            disp(['Values: ' num2str(parameter_to_vary_values_vec)])
        end
    end

    % Select whether to do one tuple of the 7 fitted parameters per country (do_geometric_median = true) or all 200 tuples of the 7 fitted parameters per country (do_geometric_median = false)
    if do_geometric_median
        assert(strcmp(stochas_run_str,'1'))
    end



    ISO_array = cell(1,num_countries); % an 1-by-n cell array of empty matrices
    country_name_array = cell(1,num_countries); % an 1-by-n cell array of empty matrices
    Regions_map = containers.Map;
    Runs_map = containers.Map;
    if ~suppress_screen_output
        if limit_to_AFRO_BD_0
            country_hours_vec = zeros(1,num_countries);
        else
            country_hours_vec = -99 * ones(1,num_countries);
        end
    end



    for country_num = 1:num_countries


        ISO = ListOfISOs{country_num};
        ISO_array{country_num} = ISO;
        demog = params_map(ISO);
        country_name_array{country_num} = demog.country_name;
        if ~suppress_screen_output
            disp([demog.country_name ' (' ISO ')']) % country's name
        end
        region_info.region_WHO = WHO_region_map(ISO);
        if ~suppress_screen_output
            disp(['WHO region: ' region_info.region_WHO])
        end
        assert(length(fields(region_info))==1)
        Regions_map(ISO) = region_info;

    
        sub_begin_time = datetime('now');
        if ~suppress_screen_output
            disp(sub_begin_time)
        end


        treatment_boundaries_vec = treatment_rates_map(ISO);
        assert(isequal(size(treatment_boundaries_vec),[1 6])) 
        % treatment year, 
        % rate to keep number of people in treatment constant, rate to have 40% of eligible people in treatment by 2030, whether the constant rate is less than the 40% rate
        % rate to have 80% of eligible people in treatment by 2030, whether the 40% rate is less than the 80% rate
        treat_start_year = treatment_boundaries_vec(1);
        assert(treat_start_year==2016)
        if do_treatment % if do_treatment is true
            in_treatment_2016_CDA = num_in_treatment_2016_map(ISO);
            assert(in_treatment_2016_CDA>=0)
            if in_treatment_2016_CDA>0
                pop_size_HBsAg_treatment_2016_vec = pop_size_HBsAg_treatment_map(ISO);
                HBsAg_treat_cov_all_ages = pop_size_HBsAg_treatment_2016_vec(4); % proportion of HBsAg+ individuals that are in treatment
                assert(HBsAg_treat_cov_all_ages>0)
                assert(HBsAg_treat_cov_all_ages<1)
                assert(in_treatment_2016_CDA==pop_size_HBsAg_treatment_2016_vec(5))
            else
                HBsAg_treat_cov_all_ages = 0;
            end
        else
            assert(strcmp(treatment_string,'no_treat'))
            HBsAg_treat_cov_all_ages = 0;
        end
    
    
        % demographic parameters:
        if ~do_geometric_median % each of the 200 stochastic particles has its own Efficacy_BirthDoseVacc_HbEAg and Efficacy_InfantVacc values, so delete the default values
            demog = rmfield(demog,'Efficacy_BirthDoseVacc_HbEAg');
            demog = rmfield(demog,'Efficacy_InfantVacc');
        end
        demog_fields = fields(demog);
        assert(ismember('Pop_byAgeGroups_1950',demog_fields))
        assert(ismember('fert',demog_fields))
        assert(ismember('MortalityRate_Women',demog_fields))
        assert(ismember('MortalityRate_Men',demog_fields))
        assert(ismember('beta_1to15',demog_fields))
        assert(ismember('beta_5plus',demog_fields))
        assert(demog.ReducInTransmission==0)
        assert(demog.YearReducInTransmission==2100)
        assert(demog.ClearanceRateWomenCoFactor==1)
        assert(ismember('CancerDeathRate',demog_fields))
        demog.PriorLAMTreatRate = 0; 
        assert(demog.Efficacy_BirthDoseVacc_HbSAg==0.95)
        assert(demog.p_VerticalTransmission_HbEAg_NoIntv==0.9)
        demog.dwvec = dwvec;

        % assign particle values
        country_start_col = country_start_cols(find(strcmp(ISO,ListOfAllISOs)));
        demog.beta_U5 = stochas_params_mat(stochas_run_num,country_start_col);
        demog.SpeedUpELoss_F = 9.5;
        demog.SpeedUpELoss_Beta = stochas_params_mat(stochas_run_num,country_start_col+1);
        demog.p_VerticalTransmission_HbSAg_NoIntv = stochas_params_mat(stochas_run_num,country_start_col+2);
        demog.cancer_rate_coeff = stochas_params_mat(stochas_run_num,country_start_col+3);
        demog.cirrh_rate_coeff = stochas_params_mat(stochas_run_num,country_start_col+4);
        demog.CancerRate_WomenCoFactor = 1;
        demog.CirrhosisRate_WomenCoFactor = 1;    
        demog.CancerRate_MenCoFactor = stochas_params_mat(stochas_run_num,country_start_col+5);
        demog.CirrhosisRate_MenCoFactor = stochas_params_mat(stochas_run_num,country_start_col+6);
        if do_treatment
            if do_geometric_median
                demog.PriorTDFTreatRate = mean(treatment_boundaries_vec(2:3));
            else
                demog.PriorTDFTreatRate = stochas_params_mat(stochas_run_num,country_start_col+7);
            end
            assert((demog.PriorTDFTreatRate>=treatment_boundaries_vec(2)) && (demog.PriorTDFTreatRate<=treatment_boundaries_vec(3)))
        else
            demog.PriorTDFTreatRate = 0;
        end
        if do_geometric_median % if doing the geometric median, use default Efficacy_BirthDoseVacc_HbEAg and Efficacy_InfantVacc values
            assert(demog.Efficacy_BirthDoseVacc_HbEAg==0.83)
            assert(demog.Efficacy_InfantVacc==0.95)
            demog.pr_VerticalTransmission_HbSAgLowVL_PAP = mean([0 0.05]); % 0.025
        else % each of the 200 stochastic particles has its own Efficacy_BirthDoseVacc_HbEAg and Efficacy_InfantVacc values
            demog.Efficacy_BirthDoseVacc_HbEAg = stochas_params_mat(stochas_run_num,end-2);
            demog.Efficacy_InfantVacc = stochas_params_mat(stochas_run_num,end-1);
            demog.pr_VerticalTransmission_HbSAgLowVL_PAP = stochas_params_mat(stochas_run_num,end);
        end



        % WUENIC coverage data released in July 2020
        infant_years = cellfun(@(yyy) str2num(erase(yyy,"x")),HepB3_table.Properties.VariableNames);
        InfantVacc_vec = HepB3_table;
        InfantVacc_vec = InfantVacc_vec(ISO,:);
        InfantVacc_vec = InfantVacc_vec{:,:};
        assert(infant_years(1)==start_year_vacc_coverage)

        BD_years = cellfun(@(yyy) str2num(erase(yyy,"x")),BD_table.Properties.VariableNames);
        BirthDose_vec = BD_table;
        BirthDose_vec = BirthDose_vec(ISO,:);
        BirthDose_vec = BirthDose_vec{:,:};
        assert(BD_years(1)==start_year_vacc_coverage)
        assert(isequal(BD_years,infant_years))
        num_years_latest_coverage = length(BD_years);
       
        demog.InfantVacc = InfantVacc_vec;
        demog.BirthDose = BirthDose_vec;
        assert(isequal(size(demog.InfantVacc),[1 (2019-1979)]))
        assert(isequal(size(demog.BirthDose),[1 (2019-1979)]))
        assert(all(demog.InfantVacc>=0))
        assert(all(demog.InfantVacc<=1))
        assert(all(demog.BirthDose>=0))
        assert(all(demog.BirthDose<=1))
    
    
        if limit_to_AFRO_BD_0
            if ~strcmp(region_info.region_WHO,'AFRO')
                continue % proceed to next country
            end
            if demog.BirthDose(end)>0
                continue % proceed to next country
            end
        end
    
    
        country_ref_data_struct = country_s_e_HCCdeaths_map(ISO);
        source_HBsAg = country_ref_data_struct.source_HBsAg;
        HBsAg_prevs_year_1 = country_ref_data_struct.HBsAg_prevs_year_1;
        demog.country_HBsAg_prevalences_by_ages_mid_1_young_old = country_ref_data_struct.country_HBsAg_prevalences_by_ages_mid_1_young_old;
        if strcmp(source_HBsAg,'Cui')
            demog.HBsAg_prevs_middle_year_1 = country_ref_data_struct.HBsAg_prevs_middle_year_1;
        else
            assert(~ismember('HBsAg_prevs_middle_year_1',fields(country_ref_data_struct)))
        end
        if strcmp(source_HBsAg,'WHO')
            demog.country_HBsAg_prevalences_by_ages_prevacc_young_old = country_ref_data_struct.country_HBsAg_prevalences_by_ages_prevacc_young_old;
        else
            assert(~ismember('country_HBsAg_prevalences_by_ages_prevacc_young_old',fields(country_ref_data_struct)))
        end



        % Progression Parameters
        % Prog gives basic (identical across age, sex and treatment) progression rates between disease states
        % Transactions then gives one finer control by allowing one to adjust rates according to age, sex and treatment
        % Transactions.Values is then used in the disease progression part of the model

        Prog = zeros(num_states, num_states); % Non-Age Specific Prog parameters stored as (from, to)

        % Fill-in transitions from Immune Tolerant
        Prog(2, 3) = 0.1;

        % Fill-in transitions from Immune Reactive
        Prog(3, 4) = 0.05;
        Prog(3, 5) = 0.005;
        Prog(3, 6) = 0.028;        % added 25.3.16

        % Fill-in constant transitions from Asymptomatic Carrier
        Prog(4, 5) = 0.01;
        Prog(4, 9) = 0.01;

        % Fill-in transitions from Chronic Hep B.
        Prog(5, 6) = 0.04;

        % Fill-in transitions from Compensated Cirrhosis
        Prog(6, 7) = 0.04;
        Prog(6, 11) = 0.04;

        % Fill-in transitions from Decompensated Cirrhsis
        Prog(7, 8) = 0.04;
        Prog(7, 11) = 0.30;

        % Fill-in transitions from Liver Cancer
        Prog(8, 11) = demog.CancerDeathRate;

        % Fill-in transition from TDF-Treatment
        Prog(10, 11) = 0.001;

        % Fill-in transition from 3TC-Treatment
        Prog(12, 13) = 0.2;

        % Fill-in transition from Failed 3TC-Treatment
        Prog(13, 8) = 0.04;
        Prog(13, 11) = 0.3;

        % Fill-in transition from Severe acute
        Prog(15, 11) = CFR_Acute * rate_6months;



        Transactions = make_transactions(Prog,demog,ages,num_age_steps,rate_6months,CFR_Acute,p_ChronicCarriage);



        effparams = make_effparams(demog);
        if ~suppress_screen_output
            disp('Finished making effparams.')
        end
        if ~isempty(parameter_to_vary)
            parameter_to_vary_original_value = effparams.(parameter_to_vary);
            assert(ismember(parameter_to_vary_original_value,parameter_to_vary_values_vec))
        end



        i84 = 85; % assummed life-expectancy at birth in all countries (all countries treated equally in terms of life-expectancy)


        % Define and run the particular scenarios:


        discount_start_year = 2024;



        Label = cell(1,num_strategies);
        Runs = cell(parameter_to_vary_num_vals,num_strategies);


        for param_val_num=1:parameter_to_vary_num_vals


            if ~isempty(parameter_to_vary)
                effparams.(parameter_to_vary) = parameter_to_vary_values_vec(param_val_num);
                % assuming that parameter_to_vary has not been used in any calculations yet!
                effparams.(other_parameter_to_vary) = parameter_to_vary_values_vec(param_val_num);
                assert(effparams.(parameter_to_vary)==effparams.(other_parameter_to_vary))
            end


            % 1: SQ - infant vaccination up to 90% and everything else at current levels continuing into the future
            strategy_num = 1;
            params = [];
            params.InfantVaccIntvCov = 0.9;
            params.TScaleup_InfantVaccineIntv_start = 2022;
            params.TScaleup_InfantVaccineIntv_finish = 2024;
    
            [Runs{param_val_num,strategy_num}] = HBVmodel(...
                demog,params,effparams,source_HBsAg,num_states,num_year_divisions,dt,ages,num_age_steps,...
                start_year_simulations,end_year_simulations,T0,start_year_vacc_coverage,first_year_results,last_year_results,num_years_results,...
                theta,ECofactor,treat_start_year-dt,HBsAg_treat_cov_all_ages,p_ChronicCarriage,Prog,Transactions);
            Runs{param_val_num,strategy_num}.DALYPerYear = make_daly_mat(Runs{param_val_num,strategy_num},num_years_results,i84);
            assert(isequal(size(Runs{param_val_num,strategy_num}.DALYPerYear),[1 num_years_results]))
            Label{strategy_num} = 'HepB3';


            Runs{param_val_num,strategy_num} = modify_run_results(i5y,Runs{param_val_num,strategy_num});


   
            if save_results


                % 2: 1 + Increase BD to 90%
                strategy_num = 2;
                if isempty(replacement_list_of_strategies) || ismember(num2str(strategy_num),replacement_list_of_strategies)
                    params = [];
                    params.InfantVaccIntvCov = 0.9;
                    params.BirthDoseIntvCov = BD_scaleup_cov;
                    params.TScaleup_InfantVaccineIntv_start = 2022;
                    params.TScaleup_InfantVaccineIntv_finish = 2024;
                    params.TScaleup_BirthDoseVaccineIntv_start = TScaleup_BirthDoseVaccineIntv_start;
                    params.TScaleup_BirthDoseVaccineIntv_finish = TScaleup_BirthDoseVaccineIntv_finish;
    
                    [Runs{param_val_num,strategy_num}] = HBVmodel(...
                        demog,params,effparams,source_HBsAg,num_states,num_year_divisions,dt,ages,num_age_steps,...
                        start_year_simulations,end_year_simulations,T0,start_year_vacc_coverage,first_year_results,last_year_results,num_years_results,...
                        theta,ECofactor,treat_start_year-dt,HBsAg_treat_cov_all_ages,p_ChronicCarriage,Prog,Transactions);
                    Runs{param_val_num,strategy_num}.DALYPerYear = make_daly_mat(Runs{param_val_num,strategy_num},num_years_results,i84);
                    Runs{param_val_num,strategy_num} = modify_run_results(i5y,Runs{param_val_num,strategy_num});

                    Label{strategy_num} = 'HepB-BD';
                end
    


                % 3: 1 + PAP to 100% of the pregnant women who have High VL (whatever eAg status) (irrespective of BD)
                strategy_num = 3;
                if isempty(replacement_list_of_strategies) || ismember(num2str(strategy_num),replacement_list_of_strategies)
                    params = [];
                    params.InfantVaccIntvCov = 0.9;
                    params.TScaleup_InfantVaccineIntv_start = 2022;
                    params.TScaleup_InfantVaccineIntv_finish = 2024;
    
                    params.cov_BirthDoseAndTDF_EAgHighVL = 1.0;
                    params.cov_BirthDoseAndTDF_SAgHighVL = 1.0;
                    params.cov_BirthDoseAndTDF_EAgLowVL = 0;      
                    params.cov_BirthDoseAndTDF_SAgLowVL = 0;
    
                    params.cov_TDFOnly_EAgHighVL=1.0;   % coverage among those who missed BD
                    params.cov_TDFOnly_SAgHighVL=1.0;
                    params.cov_TDFOnly_EAgLowVL=0;
                    params.cov_TDFOnly_SAgLowVL=0;
            
                    params.TScaleup_PAP = PAP_scaleup_end_year;
   
                    [Runs{param_val_num,strategy_num}] = HBVmodel(...
                        demog,params,effparams,source_HBsAg,num_states,num_year_divisions,dt,ages,num_age_steps,...
                        start_year_simulations,end_year_simulations,T0,start_year_vacc_coverage,first_year_results,last_year_results,num_years_results,...
                        theta,ECofactor,treat_start_year-dt,HBsAg_treat_cov_all_ages,p_ChronicCarriage,Prog,Transactions);
                    Runs{param_val_num,strategy_num}.DALYPerYear = make_daly_mat(Runs{param_val_num,strategy_num},num_years_results,i84);
                    Runs{param_val_num,strategy_num} = modify_run_results(i5y,Runs{param_val_num,strategy_num});

                    Label{strategy_num} = 'PAP-VL without HepB-BD';
                end
      


                % 4: 2 + PAP to 100% of the pregnant women who have High VL (whatever eAg status) (irrespective of BD)
                strategy_num = 4;
                if isempty(replacement_list_of_strategies) || ismember(num2str(strategy_num),replacement_list_of_strategies)
                    params = [];
                    params.InfantVaccIntvCov = 0.9;
                    params.BirthDoseIntvCov = BD_scaleup_cov;
                    params.TScaleup_InfantVaccineIntv_start = 2022;
                    params.TScaleup_InfantVaccineIntv_finish = 2024;
                    params.TScaleup_BirthDoseVaccineIntv_start = TScaleup_BirthDoseVaccineIntv_start;
                    params.TScaleup_BirthDoseVaccineIntv_finish = TScaleup_BirthDoseVaccineIntv_finish;
    
                    params.cov_BirthDoseAndTDF_EAgHighVL = 1.0;
                    params.cov_BirthDoseAndTDF_SAgHighVL = 1.0;
                    params.cov_BirthDoseAndTDF_EAgLowVL = 0;      
                    params.cov_BirthDoseAndTDF_SAgLowVL = 0;

                    params.cov_TDFOnly_EAgHighVL=1.0;   % coverage among those who missed BD
                    params.cov_TDFOnly_SAgHighVL=1.0;
                    params.cov_TDFOnly_EAgLowVL=0;
                    params.cov_TDFOnly_SAgLowVL=0;
            
                    params.TScaleup_PAP = PAP_scaleup_end_year;
   
                    [Runs{param_val_num,strategy_num}] = HBVmodel(...
                        demog,params,effparams,source_HBsAg,num_states,num_year_divisions,dt,ages,num_age_steps,...
                        start_year_simulations,end_year_simulations,T0,start_year_vacc_coverage,first_year_results,last_year_results,num_years_results,...
                        theta,ECofactor,treat_start_year-dt,HBsAg_treat_cov_all_ages,p_ChronicCarriage,Prog,Transactions);
                    Runs{param_val_num,strategy_num}.DALYPerYear = make_daly_mat(Runs{param_val_num,strategy_num},num_years_results,i84);
                    Runs{param_val_num,strategy_num} = modify_run_results(i5y,Runs{param_val_num,strategy_num});

                    Label{strategy_num} = 'PAP-VL';
                end


    
                % 5: 2 + PAP to 100% of the pregnant women who are HBeAg+ (whatever VL status) (irrespective of BD)
                strategy_num = 5;
                if isempty(replacement_list_of_strategies) || ismember(num2str(strategy_num),replacement_list_of_strategies)
                    params = [];
                    params.InfantVaccIntvCov = 0.9;
                    params.BirthDoseIntvCov = BD_scaleup_cov;
                    params.TScaleup_InfantVaccineIntv_start = 2022;
                    params.TScaleup_InfantVaccineIntv_finish = 2024;
                    params.TScaleup_BirthDoseVaccineIntv_start = TScaleup_BirthDoseVaccineIntv_start;
                    params.TScaleup_BirthDoseVaccineIntv_finish = TScaleup_BirthDoseVaccineIntv_finish;
    
                    params.cov_BirthDoseAndTDF_EAgHighVL = 1.0;
                    params.cov_BirthDoseAndTDF_SAgHighVL = 0;
                    params.cov_BirthDoseAndTDF_EAgLowVL = 1.0;      
                    params.cov_BirthDoseAndTDF_SAgLowVL = 0;    

                    params.cov_TDFOnly_EAgHighVL=1.0;   % coverage among those who missed BD
                    params.cov_TDFOnly_SAgHighVL=0;
                    params.cov_TDFOnly_EAgLowVL=1.0;
                    params.cov_TDFOnly_SAgLowVL=0;
    
                    params.TScaleup_PAP = PAP_scaleup_end_year;

                    [Runs{param_val_num,strategy_num}] = HBVmodel(...
                        demog,params,effparams,source_HBsAg,num_states,num_year_divisions,dt,ages,num_age_steps,...
                        start_year_simulations,end_year_simulations,T0,start_year_vacc_coverage,first_year_results,last_year_results,num_years_results,...
                        theta,ECofactor,treat_start_year-dt,HBsAg_treat_cov_all_ages,p_ChronicCarriage,Prog,Transactions);
                    Runs{param_val_num,strategy_num}.DALYPerYear = make_daly_mat(Runs{param_val_num,strategy_num},num_years_results,i84);
                    Runs{param_val_num,strategy_num} = modify_run_results(i5y,Runs{param_val_num,strategy_num});

                    Label{strategy_num} = 'PAP-HBeAg';    
                end




                % 6: 2 + PAP to 100% of the pregnant women who are HBsAg+, irrespective of BD
                strategy_num = 6;
                if isempty(replacement_list_of_strategies) || ismember(num2str(strategy_num),replacement_list_of_strategies)
                    params = [];
                    params.InfantVaccIntvCov = 0.9;
                    params.BirthDoseIntvCov = BD_scaleup_cov;
                    params.TScaleup_InfantVaccineIntv_start = 2022;
                    params.TScaleup_InfantVaccineIntv_finish = 2024;
                    params.TScaleup_BirthDoseVaccineIntv_start = TScaleup_BirthDoseVaccineIntv_start;
                    params.TScaleup_BirthDoseVaccineIntv_finish = TScaleup_BirthDoseVaccineIntv_finish;
    
                    params.cov_BirthDoseAndTDF_EAgHighVL = 1.0;
                    params.cov_BirthDoseAndTDF_SAgHighVL = 1.0;
                    params.cov_BirthDoseAndTDF_EAgLowVL = 1.0;      
                    params.cov_BirthDoseAndTDF_SAgLowVL = 1.0;

                    params.cov_TDFOnly_EAgHighVL=1.0;   % coverage among those who missed BD
                    params.cov_TDFOnly_SAgHighVL=1.0;
                    params.cov_TDFOnly_EAgLowVL=1.0;
                    params.cov_TDFOnly_SAgLowVL=1.0;
            
                    params.TScaleup_PAP = PAP_scaleup_end_year;
   
                    [Runs{param_val_num,strategy_num}] = HBVmodel(...
                        demog,params,effparams,source_HBsAg,num_states,num_year_divisions,dt,ages,num_age_steps,...
                        start_year_simulations,end_year_simulations,T0,start_year_vacc_coverage,first_year_results,last_year_results,num_years_results,...
                        theta,ECofactor,treat_start_year-dt,HBsAg_treat_cov_all_ages,p_ChronicCarriage,Prog,Transactions);
                    Runs{param_val_num,strategy_num}.DALYPerYear = make_daly_mat(Runs{param_val_num,strategy_num},num_years_results,i84);
                    Runs{param_val_num,strategy_num} = modify_run_results(i5y,Runs{param_val_num,strategy_num});

                    Label{strategy_num} = 'PAP-Universal';   
                end


            end % end if save_results loop


        end % end param_val_num for loop


    
        if ~suppress_screen_output
            disp('Finished running the model.')
        end



        if save_results


            assert(isequal(size(Runs),[parameter_to_vary_num_vals num_strategies]))
            assert(length(Label)==num_strategies)
            if isempty(replacement_list_of_strategies)
                assert(~(any(any(cellfun(@isempty,Runs))))) % ensure that no cell is empty
                assert(~(any(cellfun(@isempty,Label)))) % ensure that no cell is empty
            else
                strategies_all = 1:num_strategies;
                strategies_analysed = cellfun(@str2num,replacement_list_of_strategies);
                strategies_not_analysed = setdiff(strategies_all,strategies_analysed);
                assert(all(all(cellfun(@isempty,Runs(:,strategies_not_analysed))))) % ensure that all cells are empty
                assert(~(any(any(cellfun(@isempty,Runs(:,strategies_analysed)))))) % ensure that no cell is empty
                assert(all(cellfun(@isempty,Label(:,strategies_not_analysed)))) % ensure that all cells are empty
                assert(~(any(cellfun(@isempty,Label(:,strategies_analysed))))) % ensure that no cell is empty
            end
            Runs_map(ISO) = Runs;
            if do_geometric_median
                particles_string = 'geo_med';
            else
                assert(~do_geometric_median)
                particles_string = ['part_' stochas_run_str];
            end
            if do_geometric_median % only want to save Regions_map for one particle per country
                save(fullfile(folder_path_out,['Regions_map_' aaa '_countries_' particles_string '.mat']),'Regions_map')
                if ~isempty(parameter_to_vary)
                    save(fullfile(folder_path_out,'Parameter_varied_values.mat'),'parameter_to_vary','parameter_to_vary_original_value','parameter_to_vary_values_vec')
                end
            end
            save(fullfile(folder_path_out,['Runs_map_' aaa '_BD_' BD_scaleup_cov_str parameter_to_vary_short '_countries_' particles_string '.mat']),'Runs_map')
            if country_num==num_countries % labels same for all countries so only save for last country
                save(fullfile(folder_path_out,'Label_array.mat'),'Label')
            end


        end % end if save_results loop


        if ~suppress_screen_output
            time_taken_for_country = hours(datetime('now') - sub_begin_time);
            country_hours_vec(country_num) = time_taken_for_country;
            assert(all(country_hours_vec(1:country_num)>=0))
            average_time_per_country = mean(country_hours_vec(1:country_num));
            min_time_per_country = min(country_hours_vec(1:country_num));
            max_time_per_country = max(country_hours_vec(1:country_num));
            num_countries_left = num_countries - country_num;
            mean_time_left = num_countries_left * average_time_per_country;
            min_time_left = num_countries_left * min_time_per_country;
            max_time_left = num_countries_left * max_time_per_country;
            if num_countries_left>10
                if rem(num_countries_left, 10) == 0 % only show for every 10th country
                    disp([...
                        'There are ' ...
                        num2str(num_countries_left) ...
                        ' countries left, which will take about ' num2str(mean_time_left) ' (' num2str(min_time_left) ' to ' num2str(max_time_left) ') hours.'...
                        ])
                end
            else
                disp(['There are ' num2str(num_countries_left) ' countries left, which will take about ' num2str(mean_time_left) ' (' num2str(min_time_left) ' to ' num2str(max_time_left) ') hours.'])
            end
        end
    end % end for country_num loop



    if ~suppress_screen_output
        assert(length(country_hours_vec)==num_countries)

        toc

        disp('Run finished!')

        end_time = datetime('now');
        disp(end_time)
        disp(['The duration of this run was ' char(end_time - begin_time) '\n\n'])
    end


end % end function country_level_analyses
        





function effparams = make_effparams(demog)


    effparams.FracEPosHighVL = 0.90; % proportion of HBsAg+ HBeAg+ pregnant women that have high viral load
    effparams.FracSPosHighVL = 0.05; % proportion of HBsAg+ HBeAg- pregnant women that have high viral load
    p_HbSAg_av = demog.p_VerticalTransmission_HbSAg_NoIntv;
    p_HbEAg_av = 0.9; % p_VerticalTransmission_HbEAg_NoIntv


    % -------
    % ** WAY OF INTERPRETING AVAILABLE DATA #1: Transmission parameters option:
    % just differentiated by VL **
    % [Assume that transmission rate only depends on VL ..[ and make the rates 
    % by e/s-status and VL consistent with the given average of transmsision 
    % among e and s-status]

    % Given average of HBsAg is p_VerticalTransmission_HbSAg_NoIntv)
    % Given average of HBeAg is 0.90 (using in previous fitting)
    
    % get acceptable values for 'FracSPosHighVL'
    range_for_FracSPosHighVL = get_range_of_FracSPosHighVL(effparams.FracEPosHighVL, p_HbSAg_av, p_HbEAg_av);

    assert(effparams.FracSPosHighVL>=range_for_FracSPosHighVL(1))
    assert(effparams.FracSPosHighVL<=range_for_FracSPosHighVL(2))

    % Work out the transmission rates for low and high VL (for this
    % case of no differentiation by E/S state).
    % Use the constraints of the average E and S rates, and the fact that
    % transmission rate must not vary within VL categories
    beta_low_vl = max(1e-10,( p_HbEAg_av * effparams.FracSPosHighVL - p_HbSAg_av * effparams.FracEPosHighVL ) / ( effparams.FracSPosHighVL - effparams.FracEPosHighVL ) );
    beta_high_vl = (p_HbSAg_av - (1-effparams.FracSPosHighVL)*beta_low_vl)/effparams.FracSPosHighVL;
    assert(beta_low_vl>0); assert(beta_low_vl<1); assert(beta_high_vl>0); assert(beta_high_vl<1); assert(beta_high_vl>beta_low_vl)
    ratio_high_to_low = beta_high_vl  / beta_low_vl;

    effparams.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL=ratio_high_to_low;
    effparams.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL=ratio_high_to_low; 

    effparams = augment_effparams(effparams,demog);


    % CHECK THAT EVERYTHING IS CORRECT 
    % Compute p_VerticalTransmission_HbSAgHighVL_NoIntv and p_VerticalTransmission_HbSAgLowVL_NoIntv
    p_VerticalTransmission_HbSAgLowVL_NoIntv=p_HbSAg_av / ( (1-effparams.FracSPosHighVL) + effparams.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL*effparams.FracSPosHighVL ); 
    p_VerticalTransmission_HbSAgHighVL_NoIntv=min(1,effparams.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL * p_VerticalTransmission_HbSAgLowVL_NoIntv);
    
    p_VerticalTransmission_HbEAgLowVL_NoIntv=p_HbEAg_av / ( (1-effparams.FracEPosHighVL) + effparams.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL*effparams.FracEPosHighVL ); 
    p_VerticalTransmission_HbEAgHighVL_NoIntv=min(1,effparams.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL * p_VerticalTransmission_HbEAgLowVL_NoIntv);
    
    % average transmission rate among S+ correct:
    assert(abs(p_HbSAg_av - (p_VerticalTransmission_HbSAgLowVL_NoIntv*(1-effparams.FracSPosHighVL) + p_VerticalTransmission_HbSAgHighVL_NoIntv*effparams.FracSPosHighVL))<0.001)
    % p_HbSAg_av - (p_VerticalTransmission_HbSAgLowVL_NoIntv*(1-effparams.FracSPosHighVL) + p_VerticalTransmission_HbSAgHighVL_NoIntv*effparams.FracSPosHighVL)
    % = p_HbSAg_av - (p_VerticalTransmission_HbSAgLowVL_NoIntv*(1-effparams.FracSPosHighVL) + effparams.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL*p_VerticalTransmission_HbSAgLowVL_NoIntv*effparams.FracSPosHighVL)
    % = p_HbSAg_av - (p_VerticalTransmission_HbSAgLowVL_NoIntv - p_VerticalTransmission_HbSAgLowVL_NoIntv*effparams.FracSPosHighVL + effparams.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL*p_VerticalTransmission_HbSAgLowVL_NoIntv*effparams.FracSPosHighVL)
    % = p_HbSAg_av - (p_VerticalTransmission_HbSAgLowVL_NoIntv * (1 - effparams.FracSPosHighVL + effparams.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL*effparams.FracSPosHighVL))
    % = p_HbSAg_av - p_HbSAg_av
    % = 0
    
    % average transmission rate among E+ correct:
    assert(abs(p_HbEAg_av - (p_VerticalTransmission_HbEAgLowVL_NoIntv*(1-effparams.FracEPosHighVL) + p_VerticalTransmission_HbEAgHighVL_NoIntv*effparams.FracEPosHighVL))<0.001)
    % p_HbEAg_av - (p_VerticalTransmission_HbEAgLowVL_NoIntv*(1-effparams.FracEPosHighVL) + p_VerticalTransmission_HbEAgHighVL_NoIntv*effparams.FracEPosHighVL)
    % = p_HbEAg_av - (p_VerticalTransmission_HbEAgLowVL_NoIntv*(1-effparams.FracEPosHighVL) + effparams.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL*p_VerticalTransmission_HbEAgLowVL_NoIntv*effparams.FracEPosHighVL)
    % = p_HbEAg_av - (p_VerticalTransmission_HbEAgLowVL_NoIntv - p_VerticalTransmission_HbEAgLowVL_NoIntv*effparams.FracEPosHighVL + effparams.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL*p_VerticalTransmission_HbEAgLowVL_NoIntv*effparams.FracEPosHighVL)
    % = p_HbEAg_av - (p_VerticalTransmission_HbEAgLowVL_NoIntv * (1 - effparams.FracEPosHighVL + effparams.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL*effparams.FracEPosHighVL))
    % = p_HbEAg_av - p_HbEAg_av

    % equal transmssion rates within a VL category
    assert(abs(p_VerticalTransmission_HbSAgLowVL_NoIntv-p_VerticalTransmission_HbEAgLowVL_NoIntv)<0.001);

    assert(abs(p_VerticalTransmission_HbSAgHighVL_NoIntv-p_VerticalTransmission_HbEAgHighVL_NoIntv)<0.001);
    % p_VerticalTransmission_HbSAgHighVL_NoIntv-p_VerticalTransmission_HbEAgHighVL_NoIntv
    % = effparams.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL * p_VerticalTransmission_HbSAgLowVL_NoIntv - (effparams.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL * p_VerticalTransmission_HbEAgLowVL_NoIntv)
    % = ratio_high_to_low * p_VerticalTransmission_HbSAgLowVL_NoIntv - (ratio_high_to_low * p_VerticalTransmission_HbEAgLowVL_NoIntv)
    % = ratio_high_to_low * (p_VerticalTransmission_HbSAgLowVL_NoIntv - p_VerticalTransmission_HbEAgLowVL_NoIntv)
    % = ratio_high_to_low * 0 % from the previous assert
    % = 0


end % end function make_effparams






function effparams = augment_effparams(effparams,demog)

    effparams.pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc = 0.0;              % probability ratio relative to the corresponding group's _NoIntv level; 2022/03/01
    effparams.pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP = 0.0;          % probability ratio relative to the corresponding group's _NoIntv level
    effparams.pr_VerticalTransmission_HbSAgLowVL_PAP = demog.pr_VerticalTransmission_HbSAgLowVL_PAP; % probability ratio relative to the corresponding group's _NoIntv level; 2022/12/02

    effparams.pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc = 0.20;            % probability ratio relative to the corresponding group's _NoIntv level
    effparams.pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP = 0.01;        % probability ratio relative to the corresponding group's _NoIntv level
    effparams.pr_VerticalTransmission_HbSAgHighVL_PAP = 0.06;                      % probability ratio relative to the corresponding group's _NoIntv level

    effparams.pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc = effparams.pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc;
    effparams.pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP = effparams.pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP;
    effparams.pr_VerticalTransmission_HbEAgLowVL_PAP = effparams.pr_VerticalTransmission_HbSAgLowVL_PAP;

    effparams.pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc = effparams.pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc;
    effparams.pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP = effparams.pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP;
    effparams.pr_VerticalTransmission_HbEAgHighVL_PAP = effparams.pr_VerticalTransmission_HbSAgHighVL_PAP;

end % end function augment_effparams






function Transactions = make_transactions(Prog,demog,ages,num_age_steps,rate_6months,CFR_Acute,p_ChronicCarriage)

    % Summarise this matrix as transactions lists.
    [transactions_from transactions_to] = find(Prog > 0);

    % Load these into a data-structure:
    Transactions = [];
    for tr = 1:length(transactions_from)
        Transactions.From(tr) = transactions_from(tr);
        Transactions.To(tr) = transactions_to(tr);
        Transactions.Values(tr) = {repmat(Prog(transactions_from(tr), transactions_to(tr)), [1, num_age_steps, 2, 2])}; 
        % Transaction.Values is identical across age, sex, treatment
        % For each to-from pair, form a 1 x num_age_steps x 2 x 2 double containing the progression parameter for that to-from transition
        % Arrange this sequence of matrices in a cell array called Transactions.Values, which is contained in Transactions
    end


    % Add age-specific progression (14, 15)-->(2, 9)
    Transactions.From(length(Transactions.From) + 1) = 14;
    Transactions.To(length(Transactions.To) + 1) = 2;
    Transactions.Values(length(Transactions.Values) + 1) = {p_ChronicCarriage * rate_6months};

    Transactions.From(length(Transactions.From) + 1) = 14;
    Transactions.To(length(Transactions.To) + 1) = 9;
    Transactions.Values(length(Transactions.Values) + 1) = {(1 - p_ChronicCarriage) * rate_6months};

    Transactions.From(length(Transactions.From) + 1) = 15;
    Transactions.To(length(Transactions.To) + 1) = 2;
    Transactions.Values(length(Transactions.Values) + 1) = {p_ChronicCarriage * (1 - CFR_Acute) * rate_6months};

    Transactions.From(length(Transactions.From) + 1) = 15;
    Transactions.To(length(Transactions.To) + 1) = 9;
    Transactions.Values(length(Transactions.Values) + 1) = {(1 - p_ChronicCarriage) * (1 - CFR_Acute) * rate_6months};



    % Add age-specific progression (2,3) and (3,4)
    AgeSpecELossFunction = demog.SpeedUpELoss_F*exp(-demog.SpeedUpELoss_Beta*ages);

    Transactions.From(length(Transactions.From)+1) = 2;
    Transactions.To(length(Transactions.To)+1) = 3;
    Transactions.Values(length(Transactions.Values)+1) = {repmat(Prog(2,3)*AgeSpecELossFunction,[1 1 2 2])};

    Transactions.From(length(Transactions.From)+1) = 3;
    Transactions.To(length(Transactions.To)+1) = 4;
    Transactions.Values(length(Transactions.Values)+1) = {repmat(Prog(3,4)*AgeSpecELossFunction,[1 1 2 2])};



    % Modify (3,5) to an age-specific progression
    shortcut_rate_age = 20;
    indicator_vec = [zeros(1,shortcut_rate_age*10) ones(1,num_age_steps - shortcut_rate_age*10)];
    tmp_pos = find((Transactions.From == 3) & (Transactions.To == 5));
    assert(length(tmp_pos)==1)
    tmp_mat = Transactions.Values(tmp_pos);
    tmp_mat = tmp_mat{1};
    assert(min(min(min(min(tmp_mat))))==max(max(max(max(tmp_mat)))))
    Transactions.Values(tmp_pos) = {repmat(Prog(3,5)*indicator_vec,[1 1 2 2])}; % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair

    trans_rate_age = 25; % generic setting
    trans_rate_by_age = (demog.cirrh_rate_coeff*(ages - trans_rate_age)).^2; 
    trans_rate_by_age = trans_rate_by_age .* [zeros(1,trans_rate_age*10) ones(1,num_age_steps - trans_rate_age*10)];

    % Modify (5,6) to an age-specific progression
    trans_rate_by_age_5_6 = Prog(5,6)*trans_rate_by_age;
    tmp_pos = find((Transactions.From == 5) & (Transactions.To == 6));
    assert(length(tmp_pos)==1)
    tmp_mat = Transactions.Values(tmp_pos);
    tmp_mat = tmp_mat{1};
    assert(min(min(min(min(tmp_mat))))==max(max(max(max(tmp_mat)))))
    tmp_mat = repmat(trans_rate_by_age_5_6,[1 1 2 2]); % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair
    tmp_mat(:,:,1,:) = tmp_mat(:,:,1,:)*demog.CirrhosisRate_WomenCoFactor;
    tmp_mat(:,:,2,:) = tmp_mat(:,:,2,:)*demog.CirrhosisRate_MenCoFactor;
    tmp_mat = min(5, tmp_mat);
    Transactions.Values(tmp_pos) = {tmp_mat};

    % Modify (3,6) to an age-specific progression
    trans_rate_by_age_3_6 = Prog(3,6)*trans_rate_by_age;
    tmp_pos = find((Transactions.From == 3) & (Transactions.To == 6));
    assert(length(tmp_pos)==1)
    tmp_mat = Transactions.Values(tmp_pos);
    tmp_mat = tmp_mat{1};
    assert(min(min(min(min(tmp_mat))))==max(max(max(max(tmp_mat)))))
    tmp_mat = repmat(trans_rate_by_age_3_6,[1 1 2 2]); % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair
    tmp_mat(:,:,1,:) = tmp_mat(:,:,1,:)*demog.CirrhosisRate_WomenCoFactor;
    tmp_mat(:,:,2,:) = tmp_mat(:,:,2,:)*demog.CirrhosisRate_MenCoFactor;
    tmp_mat = min(5, tmp_mat);
    Transactions.Values(tmp_pos) = {tmp_mat};


    % Add sex-specific co-factor to clearance (4, 9)
    tmp_pos = find((Transactions.From == 4) & (Transactions.To == 9));
    assert(length(tmp_pos)==1)
    tmp = Transactions.Values{tmp_pos}; % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair
    tmp(:, :, 1, :) = tmp(:, :, 1, :) * demog.ClearanceRateWomenCoFactor;
    Transactions.Values(tmp_pos) = {tmp};


    % Add age-specific progresion (2, 3, 4, 5, 6)-->8
    cancer_rate_age = 10;
    cancer_rate_by_age = (demog.cancer_rate_coeff*(ages - cancer_rate_age)).^2; 
    cancer_rate_by_age = cancer_rate_by_age .* [zeros(1,cancer_rate_age*10) ones(1,num_age_steps - cancer_rate_age*10)];
    AgeSpecCancerFunction_Women = min(1, demog.CancerRate_WomenCoFactor * cancer_rate_by_age);
    AgeSpecCancerFunction_Men = min(1, demog.CancerRate_MenCoFactor * cancer_rate_by_age);

    AgeSpecificProgToCancer = zeros(1, num_age_steps, 2, 2);
    AgeSpecificProgToCancer(:, :, 1, :) = repmat(AgeSpecCancerFunction_Women, [1, 1, 1, 2]);
    AgeSpecificProgToCancer(:, :, 2, :) = repmat(AgeSpecCancerFunction_Men, [1, 1, 1, 2]);

    Transactions.From(length(Transactions.From) + 1) = 2;
    Transactions.To(length(Transactions.To) + 1) = 8;
    Transactions.Values(length(Transactions.Values) + 1) = {AgeSpecificProgToCancer};

    Transactions.From(length(Transactions.From) + 1) = 3;
    Transactions.To(length(Transactions.To) + 1) = 8;
    Transactions.Values(length(Transactions.Values) + 1) = {2 * AgeSpecificProgToCancer};

    Transactions.From(length(Transactions.From) + 1) = 4;
    Transactions.To(length(Transactions.To) + 1) = 8;
    Transactions.Values(length(Transactions.Values) + 1) = {0.5 * AgeSpecificProgToCancer};

    Transactions.From(length(Transactions.From) + 1) = 5;
    Transactions.To(length(Transactions.To) + 1) = 8;
    Transactions.Values(length(Transactions.Values) + 1) = {2 * AgeSpecificProgToCancer};

    Transactions.From(length(Transactions.From) + 1) = 6;
    Transactions.To(length(Transactions.To) + 1) = 8;
    Transactions.Values(length(Transactions.Values) + 1) = {13 * AgeSpecificProgToCancer};

end % end function make_transactions






function outmat = make_daly_mat(model_run,num_years_results,yll_max_age)

    yll_spread = squeeze(sum(sum(model_run.Prev_Deaths_1yr(:,1:(yll_max_age-1),:),1),2)); % sum over sex and age
    % prevalence of deaths because one wants to count all those that would still be alive if they had not died from HBV
    assert(isequal(size(yll_spread),[num_years_results 1]))
    

    yld_spread = squeeze(sum(sum(model_run.yld_1yr,1),2)); % years living with the disease
    % Wikipedia: YLD = I x DW x L, where I = number of incident cases in the population, DW = disability weight of specific condition, and L = average duration of the case until remission or death (years)
    % I x L = P therefore YLD = P x DW, which is done in the model
    % model_run.yld_1yr is a 2 x 100 x num_years_results matrix 
    % sum over sexes to get a 100 x 101 matrix of 100 age groups (0 to 99) versus 101 years (2000 to 2100)
    assert(isequal(size(yld_spread),[num_years_results 1]))
    outmat = yll_spread' + yld_spread';
    % DALY = YLL + YLD
    assert(isequal(size(outmat),[1 num_years_results]))

end






function run_results = modify_run_results(i5y,run_results)

    Num_preg_women_HbeAg_pos_HVL_approx = run_results.num_births_toHbEAgWomenHVL_1yr_approx; % approximated from number of babies
    Num_preg_women_HbeAg_pos_LVL_approx = run_results.num_births_toHbEAgWomenLVL_1yr_approx; % approximated from number of babies
    Num_preg_women_HbeAg_neg_HVL_approx = run_results.num_births_toHbSAgWomenHVL_1yr_approx; % approximated from number of babies
    Num_preg_women_HbeAg_neg_LVL_approx = run_results.num_births_toHbSAgWomenLVL_1yr_approx; % approximated from number of babies
    run_results.Num_preg_women_HbsAg_approx = ...
        Num_preg_women_HbeAg_pos_HVL_approx + ...
        Num_preg_women_HbeAg_pos_LVL_approx + ...
        Num_preg_women_HbeAg_neg_HVL_approx + ...
        Num_preg_women_HbeAg_neg_LVL_approx;
    run_results.Num_preg_women_HVL_approx = ...
        Num_preg_women_HbeAg_pos_HVL_approx + ...
        Num_preg_women_HbeAg_neg_HVL_approx;
    run_results.Num_preg_women_HbeAg_pos_approx = ...
        Num_preg_women_HbeAg_pos_HVL_approx + ...
        Num_preg_women_HbeAg_pos_LVL_approx;

    Num_babies_chronic_from_HbeAg_pos_HVL_preg_women = run_results.num_births_chronic_HbEAgWomenHVL_1yr_approx; % calculated from number of babies, so not approximate
    Num_babies_chronic_from_HbeAg_pos_LVL_preg_women = run_results.num_births_chronic_HbEAgWomenLVL_1yr_approx; % calculated from number of babies, so not approximate
    Num_babies_chronic_from_HbeAg_neg_HVL_preg_women = run_results.num_births_chronic_HbSAgWomenHVL_1yr_approx; % calculated from number of babies, so not approximate
    Num_babies_chronic_from_HbeAg_neg_LVL_preg_women = run_results.num_births_chronic_HbSAgWomenLVL_1yr_approx; % calculated from number of babies, so not approximate
    Num_babies_chronic = ...
        Num_babies_chronic_from_HbeAg_pos_HVL_preg_women + ...
        Num_babies_chronic_from_HbeAg_pos_LVL_preg_women + ...
        Num_babies_chronic_from_HbeAg_neg_HVL_preg_women + ...
        Num_babies_chronic_from_HbeAg_neg_LVL_preg_women;
    %run_results.Num_babies_chronic = Num_babies_chronic;
    Num_babies_HVL = ...
        Num_babies_chronic_from_HbeAg_pos_HVL_preg_women + ...
        Num_babies_chronic_from_HbeAg_neg_HVL_preg_women;
    %run_results.Num_babies_HVL = Num_babies_HVL;
    Num_babies_HbeAg = ...
        Num_babies_chronic_from_HbeAg_pos_HVL_preg_women + ...
        Num_babies_chronic_from_HbeAg_pos_LVL_preg_women;
    %run_results.Num_babies_HbeAg = Num_babies_HbeAg;

    assert(all(run_results.num_births_toHbSAgWomenHVL_1yr_approx<=run_results.num_births_toHbSAgWomenLVL_1yr_approx))
    assert(all(run_results.num_births_toHbEAgWomenLVL_1yr_approx<=run_results.num_births_toHbSAgWomenLVL_1yr_approx))
    assert(all(...
        run_results.num_births_toHbEAgWomenHVL_1yr_approx + ...
        run_results.num_births_toHbEAgWomenLVL_1yr_approx + ...
        run_results.num_births_toHbSAgWomenHVL_1yr_approx + ...
        run_results.num_births_toHbSAgWomenLVL_1yr_approx < run_results.num_births_1yr_approx))
    assert(all(run_results.num_births_chronic_HbEAgWomenHVL_1yr_approx<=run_results.num_births_toHbEAgWomenHVL_1yr_approx))
    assert(all(run_results.num_births_chronic_HbEAgWomenLVL_1yr_approx<=run_results.num_births_toHbEAgWomenLVL_1yr_approx))
    assert(all(run_results.num_births_chronic_HbSAgWomenHVL_1yr_approx<=run_results.num_births_toHbSAgWomenHVL_1yr_approx))
    assert(all(run_results.num_births_chronic_HbSAgWomenLVL_1yr_approx<=run_results.num_births_toHbSAgWomenLVL_1yr_approx))
    assert(all(run_results.PeripartumTreatment_HbEAg_HighVL_approx<=run_results.num_births_toHbEAgWomenHVL_1yr_approx))
    assert(all(run_results.PeripartumTreatment_HbEAg_LowVL_approx<=run_results.num_births_toHbEAgWomenLVL_1yr_approx))
    assert(all(run_results.PeripartumTreatment_HbSAg_HighVL_approx<=run_results.num_births_toHbSAgWomenHVL_1yr_approx))
    assert(all(run_results.PeripartumTreatment_HbSAg_LowVL_approx<=run_results.num_births_toHbSAgWomenLVL_1yr_approx))

    assert(all(abs(...
        run_results.PeripartumTreatment_HbEAg_HighVL_approx + ...
        run_results.PeripartumTreatment_HbEAg_LowVL_approx + ...
        run_results.PeripartumTreatment_HbSAg_HighVL_approx + ...
        run_results.PeripartumTreatment_HbSAg_LowVL_approx - run_results.RatePeripartumTreatment)<1e-6))


    run_results.NewChronicInfectionRate = squeeze(sum(sum(run_results.NewChronicInfectionRate,1),2))'; % 2 x 20 x num_years_results
    run_results.Tot_Pop_1yr_5_year_olds = squeeze(sum(sum(run_results.Tot_Pop_1yr(:,i5y,:),1),2))'; % 2 x 1 x num_years_results
    run_results.Tot_Pop_1yr = squeeze(sum(sum(run_results.Tot_Pop_1yr,1),2))'; % 2 x 100 x num_years_results
    run_results.NumSAg_1yr_5_year_olds = squeeze(sum(sum(run_results.NumSAg_1yr(:,i5y,:),1),2))'; % 2 x 1 x num_years_results
    %run_results.NumSAg_1yr = squeeze(sum(sum(run_results.NumSAg_1yr,1),2))'; % 2 x 100 x num_years_results
    %run_results.NumEAg_chronic_1yr = squeeze(sum(sum(run_results.NumEAg_chronic_1yr,1),2))'; % 2 x 100 x num_years_results
    %run_results.NumEAg_chronic_acute_1yr = squeeze(sum(sum(run_results.NumEAg_chronic_acute_1yr,1),2))'; % 2 x 100 x num_years_results
    run_results.Incid_Deaths_1yr_approx = squeeze(sum(sum(run_results.Incid_Deaths_1yr_approx,1),2))'; % 2 x 100 x num_years_results


    model_results_to_delete = {...
        'PeripartumTreatment_HbEAg_HighVL_approx','PeripartumTreatment_HbEAg_LowVL_approx','PeripartumTreatment_HbSAg_HighVL_approx','PeripartumTreatment_HbSAg_LowVL_approx',...
        'num_births_toHbEAgWomenHVL_1yr_approx','num_births_toHbEAgWomenLVL_1yr_approx','num_births_toHbSAgWomenHVL_1yr_approx','num_births_toHbSAgWomenLVL_1yr_approx',...
        'num_births_chronic_HbEAgWomenHVL_1yr_approx','num_births_chronic_HbEAgWomenLVL_1yr_approx','num_births_chronic_HbSAgWomenHVL_1yr_approx','num_births_chronic_HbSAgWomenLVL_1yr_approx',...
        'Tot_Pop_1yr',...
        'NumSAg_1yr','PrevEAg','NumEAg_chronic_1yr','NumEAg_chronic_acute_1yr','yld_1yr','Prev_Deaths_1yr',...
        'beta_U5','p_HbSAg_av',...
        'p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL','p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL',...
        'p_VerticalTransmission_HbSAgLowVL_NoIntv',...
        'p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc','p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP','p_VerticalTransmission_HbSAgLowVL_PAP',...
        'p_VerticalTransmission_HbSAgHighVL_NoIntv',...
        'p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc','p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP','p_VerticalTransmission_HbSAgHighVL_PAP',...
        'p_VerticalTransmission_HbEAgLowVL_NoIntv',...
        'p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc','p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP','p_VerticalTransmission_HbEAgLowVL_PAP',...
        'p_VerticalTransmission_HbEAgHighVL_NoIntv',...
        'p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc','p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP','p_VerticalTransmission_HbEAgHighVL_PAP'...
        };
    num_model_results_to_delete = length(model_results_to_delete);
    for ii=1:num_model_results_to_delete
        run_results=rmfield(run_results,model_results_to_delete{ii});
    end

    model_results_to_keep = {'Time',...
        'RateBirthDoseVacc','RateInfantVacc',...
        'PregnantWomenNeedToScreen',...
        'HBVPregnantWomenNeedToEvaluate',...
        'RatePeripartumTreatment',...
        'NumDecompCirr','NumLiverCancer',...
        'Dx_At_ANC_HBsAG','Dx_At_ANC_HBeAG','Dx_At_ANC_VL',...
        'NewChronicInfectionRate','NewChronicInfectionRate_NeonatesOnly','Tot_Pop_1yr_5_year_olds','NumSAg_1yr_5_year_olds','Incid_Deaths_1yr_approx','DALYPerYear',...
        'Num_preg_women_HbsAg_approx','Num_preg_women_HbeAg_pos_approx','Num_preg_women_HVL_approx','num_births_1yr_approx'...
        };
    assert(all(ismember(model_results_to_keep,fieldnames(run_results))))
    assert(all(ismember(fieldnames(run_results),model_results_to_keep)))


end % end function modify_run_results


