function country_level_analyses(sensitivity_analysis,...
    stochas_run_str,num_stochas_runs,...
    ListOfScenarios,num_scenarios,...
    ListOfAllISOs,ListOfISOs,num_countries,...
    BD_table,HepB3_table,...
    num_in_treatment_2016_map,pop_size_HBsAg_treatment_map,treatment_rates_map,...
    country_s_e_HCCdeaths_map,...
    params_map,dwvec,stochas_params_mat,country_start_cols,...
    filename_diaries,...
    num_states,num_year_divisions,dt,ages,num_age_steps,start_year,T0,end_year,...
    theta,CFR_Acute,rate_6months,ECofactor,p_ChronicCarriage)


    if nargin < 1
       error('No input')
    end

    assert(ismember(sensitivity_analysis,{'default','infant_100','treat_medium','treat_high'}))


    stochas_run_num = str2double(stochas_run_str);


    particles_str = ['stochastic_run_', stochas_run_str];
    filename_results = ['results_countries_', sensitivity_analysis, '_', particles_str, '.mat'];


    begin_time_run_num = datetime('now');
    disp(['Run number ' stochas_run_str ' (of ' num2str(num_stochas_runs) ') started at ' datestr(begin_time_run_num)])


    outMap = containers.Map; 
    % for this particle (stochas_run_num), outMap contains the 12 scenarios, each of which contains countryMap, which contains the model results (lastrun) for 110 countries
    label_array = cell(1,num_scenarios);
    scenario_hours_vec = repmat(duration(0,0,0),1,num_scenarios);

    for scenario_num = 1:num_scenarios

        begin_time_scenario = datetime('now');

        scenario = ListOfScenarios{scenario_num};
        %disp(scenario)
        disp(['The "' scenario '" scenario (run number ' stochas_run_str ') started at ' datestr(begin_time_scenario)])
                
        diary off
        diary(filename_diaries)


        countryMap = containers.Map; 
        % for this particle (stochas_run_num) and scenario (scenario), countryMap contains the model results (lastrun) of each of the 110 countries in this scenario


        for country_num = 1:num_countries

            ISO = ListOfISOs{country_num};

            treatment_boundaries_vec = treatment_rates_map(ISO);
            assert(isequal(size(treatment_boundaries_vec),[1 6])) 
            % treatment year, 
            % rate to keep number of people in treatment constant, rate to have 40% of eligible people in treatment by 2030, whether the constant rate is less than the 40% rate
            % rate to have 80% of eligible people in treatment by 2030, whether the 40% rate is less than the 80% rate
            treat_start_year = treatment_boundaries_vec(1);
            assert(treat_start_year==2016)
            in_treatment_2016_CDA = num_in_treatment_2016_map(ISO);
            assert(in_treatment_2016_CDA>=0)
            if in_treatment_2016_CDA>0
                pop_size_HBsAg_treatment_2016_vec = pop_size_HBsAg_treatment_map(ISO);
                HBsAg_treat_cov_all_ages = pop_size_HBsAg_treatment_2016_vec(4);
                assert(HBsAg_treat_cov_all_ages>0)
                assert(HBsAg_treat_cov_all_ages<1)
                assert(in_treatment_2016_CDA==pop_size_HBsAg_treatment_2016_vec(5))
            else
                HBsAg_treat_cov_all_ages = 0;
            end

            params = params_map(ISO);
            params = rmfield(params,'Efficacy_BirthDoseVacc_HbEAg');
            params = rmfield(params,'Efficacy_InfantVacc');

            % Parameters
            assert(params.Efficacy_BirthDoseVacc_HbSAg==0.95)
            assert(params.p_VerticalTransmission_HbEAg_NoIntv==0.9)
            params.dwvec = dwvec;

            country_start_col = country_start_cols(find(strcmp(ISO,ListOfAllISOs)));
            params.beta_U5 = stochas_params_mat(stochas_run_num,country_start_col);
            params.SpeedUpELoss_F = 9.5;
            params.SpeedUpELoss_Beta = stochas_params_mat(stochas_run_num,country_start_col+1);
            params.p_VerticalTransmission_HbSAg_NoIntv = stochas_params_mat(stochas_run_num,country_start_col+2);
            params.cancer_rate_coeff = stochas_params_mat(stochas_run_num,country_start_col+3);
            params.cirrh_rate_coeff = stochas_params_mat(stochas_run_num,country_start_col+4);
            params.CancerRate_WomenCoFactor = 1;
            params.CirrhosisRate_WomenCoFactor = 1;    
            params.CancerRate_MenCoFactor = stochas_params_mat(stochas_run_num,country_start_col+5);
            params.CirrhosisRate_MenCoFactor = stochas_params_mat(stochas_run_num,country_start_col+6);
            switch sensitivity_analysis
                case {'default','infant_100'}
                    params.PriorTDFTreatRate = stochas_params_mat(stochas_run_num,country_start_col+7);
                    assert((params.PriorTDFTreatRate>=treatment_boundaries_vec(2)) && (params.PriorTDFTreatRate<=treatment_boundaries_vec(3)))
                case 'treat_medium'
                    params.PriorTDFTreatRate = treatment_boundaries_vec(3); % 40%
                case 'treat_high'
                    params.PriorTDFTreatRate = treatment_boundaries_vec(5); % 80%
            end
            params.Efficacy_BirthDoseVacc_HbEAg = stochas_params_mat(stochas_run_num,end-1);
            params.Efficacy_InfantVacc = stochas_params_mat(stochas_run_num,end);

    

            % WUENIC coverage data released in July 2020

            infant_years = cellfun(@(yyy) str2num(erase(yyy,"x")),HepB3_table.Properties.VariableNames);
            InfantVacc_vec = HepB3_table;
            InfantVacc_vec = InfantVacc_vec(ISO,:);
            InfantVacc_vec = InfantVacc_vec{:,:};
            assert(infant_years(1)==1980)

            BD_years = cellfun(@(yyy) str2num(erase(yyy,"x")),BD_table.Properties.VariableNames);
            BirthDose_vec = BD_table;
            BirthDose_vec = BirthDose_vec(ISO,:);
            BirthDose_vec = BirthDose_vec{:,:};
            assert(BD_years(1)==1980)
            assert(isequal(BD_years,infant_years))
       
            InfantVacc = InfantVacc_vec;
            BirthDose = BirthDose_vec;
            assert(isequal(size(InfantVacc),[1 (2019-1979)]))
            assert(isequal(size(BirthDose),[1 (2019-1979)]))
            assert(all(InfantVacc>=0))
            assert(all(InfantVacc<=1))
            assert(all(BirthDose>=0))
            assert(all(BirthDose<=1))

            years_vec_01yr = start_year:dt:end_year;
            last_available_year = 2019.0;
            index_last_available_year = last_available_year - 1979;
            last_BD_scaleup_year = 2030.0;
            cov_InfantVacc = InfantVacc(1:index_last_available_year); % coverage of vaccination from 1980 to last_available_year
            cov_BirthDose = BirthDose(1:index_last_available_year); 
            % 2019 is the last year of vaccination available from WUENIC
            first_expansion_year = 2020.0;


            max_BD_cov_val_25 = max([cov_BirthDose(end) 0.25]);
            max_BD_cov_val_50 = max([cov_BirthDose(end) 0.5]);
            max_BD_cov_val_75 = max([cov_BirthDose(end) 0.75]);
            max_BD_cov_val_90 = max([cov_BirthDose(end) 0.9]);


            if strcmp(sensitivity_analysis,'infant_100')
                future_xvals_vec = [2019.0 first_expansion_year (first_expansion_year+0.1) end_year];
                future_yvals_vec = [cov_InfantVacc(end) cov_InfantVacc(end) 1 1];
                infant_vacc_vec = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_InfantVacc,future_xvals_vec,future_yvals_vec);
                infant_vacc_vec = min(1,infant_vacc_vec);
            else
                future_xvals_vec = [2019.0 first_expansion_year end_year];
                future_yvals_vec = [cov_InfantVacc(end) cov_InfantVacc(end) cov_InfantVacc(end)];
                infant_vacc_vec = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_InfantVacc,future_xvals_vec,future_yvals_vec);
                infant_vacc_vec = min(1,infant_vacc_vec);
            end
            assert(isequal(size(infant_vacc_vec),size(years_vec_01yr)))


            switch scenario
                case 'Status quo infant & BD'
                    % Current level for each country is maintained forever

                    params.InfantVacc = infant_vacc_vec;

                    future_xvals_vec = [2019.0 first_expansion_year end_year];
                    future_yvals_vec = [cov_BirthDose(end) cov_BirthDose(end) cov_BirthDose(end)];
                    params.BirthDose = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_BirthDose,future_xvals_vec,future_yvals_vec);
                    params.BirthDose = min(1,params.BirthDose);
                    assert(isequal(size(params.BirthDose),size(years_vec_01yr)))

                    label_array{scenario_num} = 'Status quo HepB3 & HepB-BD (baseline)';
                case 'Status quo infant & BD expansion to 25%'
                    % Any country that is <25%, goes to 25% by last_BD_scaleup_year (linear expansion); others that are >25% stay at that level.

                    params.InfantVacc = infant_vacc_vec;

                    future_xvals_vec = [2019.0 first_expansion_year last_BD_scaleup_year end_year];
                    future_yvals_vec = [cov_BirthDose(end) cov_BirthDose(end) max_BD_cov_val_25 max_BD_cov_val_25];
                    params.BirthDose = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_BirthDose,future_xvals_vec,future_yvals_vec);
                    params.BirthDose = min(1,params.BirthDose);
                    assert(isequal(size(params.BirthDose),size(years_vec_01yr)))

                    label_array{scenario_num} = 'HepB-BD scale-up to 25%';
                case 'Status quo infant & BD expansion to 50%'
                    % Any country that is <50%, goes to 50% by last_BD_scaleup_year (linear expansion); others that are >50% stay at that level.

                    params.InfantVacc = infant_vacc_vec;

                    future_xvals_vec = [2019.0 first_expansion_year last_BD_scaleup_year end_year];
                    future_yvals_vec = [cov_BirthDose(end) cov_BirthDose(end) max_BD_cov_val_50 max_BD_cov_val_50];
                    params.BirthDose = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_BirthDose,future_xvals_vec,future_yvals_vec);
                    params.BirthDose = min(1,params.BirthDose);
                    assert(isequal(size(params.BirthDose),size(years_vec_01yr)))

                    label_array{scenario_num} = 'HepB-BD scale-up to 50%';
                case 'Status quo infant & BD expansion to 75%'
                    % Any country that is <75%, goes to 75% by last_BD_scaleup_year (linear expansion); others that are >75% stay at that level.

                    params.InfantVacc = infant_vacc_vec;

                    future_xvals_vec = [2019.0 first_expansion_year last_BD_scaleup_year end_year];
                    future_yvals_vec = [cov_BirthDose(end) cov_BirthDose(end) max_BD_cov_val_75 max_BD_cov_val_75];
                    params.BirthDose = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_BirthDose,future_xvals_vec,future_yvals_vec);
                    params.BirthDose = min(1,params.BirthDose);
                    assert(isequal(size(params.BirthDose),size(years_vec_01yr)))

                    label_array{scenario_num} = 'HepB-BD scale-up to 75%';
                case 'Status quo infant & BD expansion to 90%'
                    % Any country that is <90%, goes to 90% by last_BD_scaleup_year (linear expansion); others that are >90% stay at that level.

                    params.InfantVacc = infant_vacc_vec;

                    future_xvals_vec = [2019.0 first_expansion_year last_BD_scaleup_year end_year];
                    future_yvals_vec = [cov_BirthDose(end) cov_BirthDose(end) max_BD_cov_val_90 max_BD_cov_val_90];
                    params.BirthDose = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_BirthDose,future_xvals_vec,future_yvals_vec);
                    params.BirthDose = min(1,params.BirthDose);
                    assert(isequal(size(params.BirthDose),size(years_vec_01yr)))

                    label_array{scenario_num} = 'HepB-BD scale-up to 90%';
                case 'Status quo infant & BD drop 5 2020'
                    % Follows status quo but, during the year 2020, birth dose vaccination drops by 5% in relative terms.

                    params.InfantVacc = infant_vacc_vec;

                    future_xvals_vec = [2019.0 2019.9 2020.0 2020.9 2021.0 2030 end_year];
                    future_yvals_vec = [cov_BirthDose(end) cov_BirthDose(end) 0.95*cov_BirthDose(end) 0.95*cov_BirthDose(end) cov_BirthDose(end) cov_BirthDose(end) cov_BirthDose(end)];
                    params.BirthDose = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_BirthDose,future_xvals_vec,future_yvals_vec);
                    params.BirthDose = min(1,params.BirthDose);
                    assert(isequal(size(params.BirthDose),size(years_vec_01yr)))
 
                    label_array{scenario_num} = 'HepB-BD disruptions 5% in 2020';
                case 'Status quo infant & BD drop 10 2020'
                    % Follows status quo but, during the year 2020, birth dose vaccination drops by 10% in relative terms.

                    params.InfantVacc = infant_vacc_vec;

                    future_xvals_vec = [2019.0 2019.9 2020.0 2020.9 2021.0 2030 end_year];
                    future_yvals_vec = [cov_BirthDose(end) cov_BirthDose(end) 0.9*cov_BirthDose(end) 0.9*cov_BirthDose(end) cov_BirthDose(end) cov_BirthDose(end) cov_BirthDose(end)];
                    params.BirthDose = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_BirthDose,future_xvals_vec,future_yvals_vec);
                    params.BirthDose = min(1,params.BirthDose);
                    assert(isequal(size(params.BirthDose),size(years_vec_01yr)))

                    label_array{scenario_num} = 'HepB-BD disruptions 10% in 2020';
                case 'Status quo infant & BD drop 15 2020'
                    % Follows status quo but, during the year 2020, birth dose vaccination drops by 15% in relative terms.

                    params.InfantVacc = infant_vacc_vec;

                    future_xvals_vec = [2019.0 2019.9 2020.0 2020.9 2021.0 2030 end_year];
                    future_yvals_vec = [cov_BirthDose(end) cov_BirthDose(end) 0.85*cov_BirthDose(end) 0.85*cov_BirthDose(end) cov_BirthDose(end) cov_BirthDose(end) cov_BirthDose(end)];
                    params.BirthDose = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_BirthDose,future_xvals_vec,future_yvals_vec);
                    params.BirthDose = min(1,params.BirthDose);
                    assert(isequal(size(params.BirthDose),size(years_vec_01yr)))

                    label_array{scenario_num} = 'HepB-BD disruptions 15% in 2020';
                case 'Status quo infant & BD drop 20 2020'
                    % Follows status quo but, during the year 2020, birth dose vaccination drops by 20% in relative terms.

                    params.InfantVacc = infant_vacc_vec;

                    future_xvals_vec = [2019.0 2019.9 2020.0 2020.9 2021.0 2030 end_year];
                    future_yvals_vec = [cov_BirthDose(end) cov_BirthDose(end) 0.8*cov_BirthDose(end) 0.8*cov_BirthDose(end) cov_BirthDose(end) cov_BirthDose(end) cov_BirthDose(end)];
                    params.BirthDose = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_BirthDose,future_xvals_vec,future_yvals_vec);
                    params.BirthDose = min(1,params.BirthDose);
                    assert(isequal(size(params.BirthDose),size(years_vec_01yr)))

                    label_array{scenario_num} = 'HepB-BD disruptions 20% in 2020';
                case 'Status quo infant & BD delayed expansion 2023 to 2030'
                    % Planned expansion of birth-dose vaccination is postponed by 3 years

                    params.InfantVacc = infant_vacc_vec;

                    future_xvals_vec = [2019.0 first_expansion_year first_expansion_year+3 last_BD_scaleup_year end_year];
                    future_yvals_vec = [cov_BirthDose(end) cov_BirthDose(end) cov_BirthDose(end) max_BD_cov_val_90 max_BD_cov_val_90];
                    params.BirthDose = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_BirthDose,future_xvals_vec,future_yvals_vec);
                    params.BirthDose = min(1,params.BirthDose);
                    assert(isequal(size(params.BirthDose),size(years_vec_01yr)))

                    label_array{scenario_num} = 'HepB-BD delayed & fast scale-up 2023 to 2030';
                case 'Status quo infant & BD delayed expansion 2023 to 2033'
                    % Planned expansion of birth-dose vaccination is postponed by 3 years

                    params.InfantVacc = infant_vacc_vec;

                    future_xvals_vec = [2019.0 first_expansion_year first_expansion_year+3 last_BD_scaleup_year+3 end_year];
                    future_yvals_vec = [cov_BirthDose(end) cov_BirthDose(end) cov_BirthDose(end) max_BD_cov_val_90 max_BD_cov_val_90];
                    params.BirthDose = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_BirthDose,future_xvals_vec,future_yvals_vec);
                    params.BirthDose = min(1,params.BirthDose);
                    assert(isequal(size(params.BirthDose),size(years_vec_01yr)))

                    label_array{scenario_num} = 'HepB-BD delayed & normal scale-up 2023 to 2033';
                case 'Status quo infant & BD delayed expansion 2025 to 2040'
                    % Planned expansion of birth-dose vaccination is postponed by 5 years and takes longer

                    params.InfantVacc = infant_vacc_vec;

                    future_xvals_vec = [2019.0 first_expansion_year first_expansion_year+5 last_BD_scaleup_year+10 end_year];
                    future_yvals_vec = [cov_BirthDose(end) cov_BirthDose(end) cov_BirthDose(end) max_BD_cov_val_90 max_BD_cov_val_90];
                    params.BirthDose = make_coverage_vec(start_year,num_year_divisions,dt,end_year,cov_BirthDose,future_xvals_vec,future_yvals_vec);
                    params.BirthDose = min(1,params.BirthDose);
                    assert(isequal(size(params.BirthDose),size(years_vec_01yr)))

                    label_array{scenario_num} = 'HepB-BD delayed & slow scale-up 2025 to 2040';
            end

    
            country_ref_data_struct = country_s_e_HCCdeaths_map(ISO);
            source_HBsAg = country_ref_data_struct.source_HBsAg;
            HBsAg_prevs_year_1 = country_ref_data_struct.HBsAg_prevs_year_1;
            params.country_HBsAg_prevalences_by_ages_mid_1_young_old = country_ref_data_struct.country_HBsAg_prevalences_by_ages_mid_1_young_old;
            if strcmp(source_HBsAg,'Cui')
                params.HBsAg_prevs_middle_year_1 = country_ref_data_struct.HBsAg_prevs_middle_year_1;
            else
                assert(~ismember('HBsAg_prevs_middle_year_1',fields(country_ref_data_struct)))
            end
            if strcmp(source_HBsAg,'WHO')
                params.country_HBsAg_prevalences_by_ages_prevacc_young_old = country_ref_data_struct.country_HBsAg_prevalences_by_ages_prevacc_young_old;
            else
                assert(~ismember('country_HBsAg_prevalences_by_ages_prevacc_young_old',fields(country_ref_data_struct)))
            end



            % Progression Parameters
            % Prog gives basic (identical across age, gender and treatment) progression rates between disease states
            % Transactions then gives one finer control by allowing one to adjust rates
            % according to age, gender and treatment
            % Transactions.Values is then used in the disease progression part of the
            % model

            Prog = zeros(num_states, num_states); % Non-Age Specific Prog parameters stored as (from, to)

            % Fill-in transitions from Immune Tolerant
            Prog(2, 3) = 0.1;

            % Fill-in transitions from Immune Reactive
            Prog(3, 4) = 0.05;
            Prog(3, 5) = 0.005;
            Prog(3, 6) = 0.028;

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
            Prog(8, 11) = params.CancerDeathRate;

            % Fill-in transition from TDF-Treatment
            Prog(10, 11) = 0.001;

            % Fill-in transition from 3TC-Treatment
            Prog(12, 13) = 0.2;

            % Fill-in transition from Failed 3TC-Treatment
            Prog(13, 8) = 0.04;
            Prog(13, 11) = 0.3;

            % Fill-in transition from Severe acute
            Prog(15, 11) = CFR_Acute * rate_6months;



            % Summarise this matrix as transactions lists.
            [transactions_from transactions_to] = find(Prog > 0);

            % Load these into a data-structure:
            Transactions = [];
            for tr = 1:length(transactions_from)
                Transactions.From(tr) = transactions_from(tr);
                Transactions.To(tr) = transactions_to(tr);
                Transactions.Values(tr) = {repmat(Prog(transactions_from(tr), transactions_to(tr)), [1, num_age_steps, 2, 2])}; 
                % Transaction.Values is identical across age, gender, treatment
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
            AgeSpecELossFunction=params.SpeedUpELoss_F*exp(-params.SpeedUpELoss_Beta*ages);

            Transactions.From(length(Transactions.From)+1)=2;
            Transactions.To(length(Transactions.To)+1)=3;
            Transactions.Values(length(Transactions.Values)+1)={repmat(Prog(2,3)*AgeSpecELossFunction,[1 1 2 2])};

            Transactions.From(length(Transactions.From)+1)=3;
            Transactions.To(length(Transactions.To)+1)=4;
            Transactions.Values(length(Transactions.Values)+1)={repmat(Prog(3,4)*AgeSpecELossFunction,[1 1 2 2])};



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
            trans_rate_by_age = (params.cirrh_rate_coeff*(ages - trans_rate_age)).^2; 
            trans_rate_by_age = trans_rate_by_age .* [zeros(1,trans_rate_age*10) ones(1,num_age_steps - trans_rate_age*10)];

            % Modify (5,6) to an age-specific progression
            trans_rate_by_age_5_6 = Prog(5,6)*trans_rate_by_age;
            tmp_pos = find((Transactions.From == 5) & (Transactions.To == 6));
            assert(length(tmp_pos)==1)
            tmp_mat = Transactions.Values(tmp_pos);
            tmp_mat = tmp_mat{1};
            assert(min(min(min(min(tmp_mat))))==max(max(max(max(tmp_mat)))))
            tmp_mat = repmat(trans_rate_by_age_5_6,[1 1 2 2]); % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair
            tmp_mat(:,:,1,:) = tmp_mat(:,:,1,:)*params.CirrhosisRate_WomenCoFactor;
            tmp_mat(:,:,2,:) = tmp_mat(:,:,2,:)*params.CirrhosisRate_MenCoFactor;
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
            tmp_mat(:,:,1,:) = tmp_mat(:,:,1,:)*params.CirrhosisRate_WomenCoFactor;
            tmp_mat(:,:,2,:) = tmp_mat(:,:,2,:)*params.CirrhosisRate_MenCoFactor;
            tmp_mat = min(5, tmp_mat);
            Transactions.Values(tmp_pos) = {tmp_mat};


            % Add sex-specific co-factor to clearance (4, 9)
            tmp_pos = find((Transactions.From == 4) & (Transactions.To == 9));
            assert(length(tmp_pos)==1)
            tmp = Transactions.Values{tmp_pos}; % 1 x num_age_steps x 2 x 2 double giving progression rates for this to-from pair
            tmp(:, :, 1, :) = tmp(:, :, 1, :) * params.ClearanceRateWomenCoFactor;
            Transactions.Values(tmp_pos) = {tmp};


            % Add age-specific progresion (2, 3, 4, 5, 6)-->8
            cancer_rate_age = 10;
            cancer_rate_by_age = (params.cancer_rate_coeff*(ages - cancer_rate_age)).^2; 
            cancer_rate_by_age = cancer_rate_by_age .* [zeros(1,cancer_rate_age*10) ones(1,num_age_steps - cancer_rate_age*10)];
            AgeSpecCancerFunction_Women = min(1, params.CancerRate_WomenCoFactor * cancer_rate_by_age);
            AgeSpecCancerFunction_Men = min(1, params.CancerRate_MenCoFactor * cancer_rate_by_age);

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



            i84 = 85; % the approximate life-expectancy at birth for persons in Japan and the highest life-expectnacy in the world, which is often taken as the benchmark


            num_year_1980_2100 = 2100 - 1980 + 1;


            %% Run scenarios:
            lastrun = HBVmodel(source_HBsAg,...
                num_states,num_year_divisions,dt,ages,num_age_steps,start_year,T0,...
                theta,ECofactor,treat_start_year-dt,HBsAg_treat_cov_all_ages,params,p_ChronicCarriage,Prog,Transactions);
            lastrun.DALYPerYear = make_daly_mat(lastrun,T0,num_year_1980_2100,i84);
            assert(isequal(size(lastrun.DALYPerYear),[100 T0 + 1]))

            assert(isequal(size(lastrun.Time),[1 (T0 + 1)]))
            i1980 = find(lastrun.Time>=1980, 1);
            i2100 = find(lastrun.Time>=2100, 1);
            num_cols_out = i2100 - i1980 + 1; % every output should have entries for the years 1980 to 2100
            i5y = 6;
            index_under_5y = 1:5;
            num_cols_in = i2100 - i1980 + 1;
            index_last_year = i2100;
            
            lastrun = rmfield(lastrun,'Prev_Deaths_1yr');
            lastrun = rmfield(lastrun,'yld_1yr');
            lastrun_fields = sort(fields(lastrun));
            fields_of_interest = {'Time',...
                'Tot_Pop_1yr','num_births_1yr',...
                'Incid_chronic_all_1yr_approx',...
                'NumSAg_1yr','NumSAg_chronic_1yr',...
                'Prev_TDF_treat_1yr','Prev_Immune_Reactive_1yr','Prev_Chronic_Hep_B_1yr','Prev_Comp_Cirr_1yr','Prev_Decomp_Cirr_1yr',...
                'Incid_Deaths_1yr_approx',...
                'DALYPerYear'};
            assert(all(ismember(fields_of_interest,lastrun_fields)))
            assert(all(ismember(lastrun_fields,fields_of_interest)))

            lastrun.Prev_treatment_eligible_1yr = ...
                lastrun.Prev_Immune_Reactive_1yr + lastrun.Prev_Chronic_Hep_B_1yr + lastrun.Prev_Comp_Cirr_1yr + lastrun.Prev_Decomp_Cirr_1yr + ...
                lastrun.Prev_TDF_treat_1yr;
            lastrun = rmfield(lastrun,'Prev_Immune_Reactive_1yr');
            lastrun = rmfield(lastrun,'Prev_Chronic_Hep_B_1yr');
            lastrun = rmfield(lastrun,'Prev_Comp_Cirr_1yr');
            lastrun = rmfield(lastrun,'Prev_Decomp_Cirr_1yr');

            assert(isequal(size(lastrun.Incid_Deaths_1yr_approx),[2 100 (T0 + 1)]))
            birth_cohorts_fun = @(yy) arrayfun(@(xx) sum(diag(yy,xx)), 0:(num_cols_out-1));
            assert(isequal(size(squeeze(sum(lastrun.Incid_Deaths_1yr_approx(:,:,i1980:index_last_year),1))),[100 num_cols_in]))
            lastrun.Incid_Deaths_1yr_approx_birth_cohorts = birth_cohorts_fun(squeeze(sum(lastrun.Incid_Deaths_1yr_approx(:,:,i1980:index_last_year),1))); % 2 x 100 x (T0 + 1)
            assert(isequal(size(lastrun.Incid_Deaths_1yr_approx_birth_cohorts),[1 num_cols_out]))
            lastrun.Time = lastrun.Time(i1980:i2100); % 1 x (T0 + 1)
            lastrun.Tot_Pop_1yr_5_year_olds = squeeze(sum(sum(lastrun.Tot_Pop_1yr(:,i5y,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)
            lastrun.Tot_Pop_1yr_under_5_year_olds = squeeze(sum(sum(lastrun.Tot_Pop_1yr(:,index_under_5y,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)
            lastrun.NumSAg_1yr_5_year_olds = squeeze(sum(sum(lastrun.NumSAg_1yr(:,i5y,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)
            lastrun.NumSAg_1yr_under_5_year_olds = squeeze(sum(sum(lastrun.NumSAg_1yr(:,index_under_5y,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)
            lastrun.NumSAg_chronic_1yr_5_year_olds = squeeze(sum(sum(lastrun.NumSAg_chronic_1yr(:,i5y,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)            
            lastrun.NumSAg_chronic_1yr_under_5_year_olds = squeeze(sum(sum(lastrun.NumSAg_chronic_1yr(:,index_under_5y,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)            
            lastrun.Tot_Pop_1yr = squeeze(sum(sum(lastrun.Tot_Pop_1yr(:,:,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)
            lastrun.num_births_1yr = lastrun.num_births_1yr(i1980:i2100); % 1 x (T0 + 1)
            lastrun.Incid_chronic_all_1yr_approx = squeeze(sum(sum(lastrun.Incid_chronic_all_1yr_approx(:,:,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)
            lastrun.NumSAg_1yr = squeeze(sum(sum(lastrun.NumSAg_1yr(:,:,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)
            lastrun.NumSAg_chronic_1yr = squeeze(sum(sum(lastrun.NumSAg_chronic_1yr(:,:,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)
            lastrun.Prev_TDF_treat_1yr = squeeze(sum(sum(lastrun.Prev_TDF_treat_1yr(:,:,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)
            lastrun.Prev_treatment_eligible_1yr = squeeze(sum(sum(lastrun.Prev_treatment_eligible_1yr(:,:,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)
            lastrun.Incid_Deaths_1yr_approx = squeeze(sum(sum(lastrun.Incid_Deaths_1yr_approx(:,:,i1980:i2100),1),2))'; % 2 x 100 x (T0 + 1)
            lastrun.DALYPerYear = squeeze(sum(lastrun.DALYPerYear(:,i1980:i2100),1)); % 100 x (T0 + 1)
            assert(isequal(size(lastrun.Time),[1 num_cols_out]))
            assert(isequal(size(lastrun.Tot_Pop_1yr_5_year_olds),[1 num_cols_out]))
            assert(isequal(size(lastrun.Tot_Pop_1yr_under_5_year_olds),[1 num_cols_out]))
            assert(isequal(size(lastrun.Tot_Pop_1yr),[1 num_cols_out]))
            assert(isequal(size(lastrun.DALYPerYear),[1 num_cols_out]))
            assert(all(lastrun.Tot_Pop_1yr>=lastrun.NumSAg_1yr))
            assert(all(lastrun.NumSAg_1yr>=lastrun.NumSAg_chronic_1yr))
            assert(all(lastrun.NumSAg_chronic_1yr>=lastrun.Prev_treatment_eligible_1yr))
            assert(all(lastrun.Prev_treatment_eligible_1yr>=lastrun.Prev_TDF_treat_1yr))
            assert(all(lastrun.Tot_Pop_1yr>=lastrun.Incid_Deaths_1yr_approx))
            assert(all(lastrun.Tot_Pop_1yr>=lastrun.Tot_Pop_1yr_5_year_olds))
            assert(all(lastrun.Tot_Pop_1yr>=lastrun.Tot_Pop_1yr_under_5_year_olds))
            assert(all(lastrun.Tot_Pop_1yr_5_year_olds>=lastrun.NumSAg_1yr_5_year_olds))
            assert(all(lastrun.Tot_Pop_1yr_under_5_year_olds>=lastrun.NumSAg_1yr_under_5_year_olds))
            assert(all(lastrun.NumSAg_1yr_5_year_olds>=lastrun.NumSAg_chronic_1yr_5_year_olds))
            assert(all(lastrun.NumSAg_1yr_under_5_year_olds>=lastrun.NumSAg_chronic_1yr_under_5_year_olds))


            if strcmp(stochas_run_str,'1')
                lastrun.country_name = params.country_name;
                lastrun.scenario = scenario;

                i1980 = find(years_vec_01yr>=1980,1);
                i2100 = find(years_vec_01yr>=2100,1);
                lastrun.InfantVacc = params.InfantVacc(i1980:i2100);
                lastrun.BirthDoseVacc = params.BirthDose(i1980:i2100);
            end


            countryMap(ISO) = lastrun;

        end % end for country_num loop

        outMap(scenario) = countryMap;
        save(filename_results,'outMap') 
        if strcmp(stochas_run_str,'1') && strcmp(sensitivity_analysis,'default')
            save('scenarios_array.mat','label_array') % only version saved after last scenario is correct
        end

        time_taken_for_scenario = datetime('now') - begin_time_scenario;
        scenario_hours_vec(scenario_num) = time_taken_for_scenario;
        assert(all(scenario_hours_vec(1:scenario_num)>0))
        average_time_per_scenario = mean(scenario_hours_vec(1:scenario_num));
        min_time_per_scenario = min(scenario_hours_vec(1:scenario_num));
        max_time_per_scenario = max(scenario_hours_vec(1:scenario_num));
        num_scenarios_left = num_scenarios - scenario_num;
        mean_time_left = num_scenarios_left * average_time_per_scenario;
        min_time_left = num_scenarios_left * min_time_per_scenario;
        max_time_left = num_scenarios_left * max_time_per_scenario;
        if num_scenarios_left>0
            disp(['There are ' num2str(num_scenarios_left) ' scenarios left for run number ' stochas_run_str ' (' sensitivity_analysis '), which will take about ' char(mean_time_left) ' hours.'])
        end


    end % end for scenario_num loop

    assert(length(scenario_hours_vec)==num_scenarios)
    end_time_run_num = datetime('now');
    disp(end_time_run_num)
    time_taken_for_run = end_time_run_num - begin_time_run_num;
    disp(['The duration of run number ' stochas_run_str ' (of ' num2str(num_stochas_runs) ') was ' char(time_taken_for_run) ' hours.\n\n'])
    if stochas_run_num<num_stochas_runs
        num_runs_left = num_stochas_runs - stochas_run_num;
        approximate_time_left = num_runs_left * time_taken_for_run;
        disp(['There are ' num2str(num_runs_left) ' runs left for sensitivity analysis "' sensitivity_analysis '", which will take about ' char(approximate_time_left) ' hours.'])
    end

end % end function country_level_analyses




function outmat = make_daly_mat(model_run,T0,num_year_1980_2100,yll_max_age)

    yll_spread = squeeze(sum(model_run.Prev_Deaths_1yr(:,1:(yll_max_age-1),:),1)); % sum over gender
    % prevalence of deaths because one wants to count total deaths
    assert(isequal(size(yll_spread),[(yll_max_age-1) T0+1]))
    yll_spread = [yll_spread; zeros(100-(yll_max_age-1),T0+1)];
    assert(isequal(size(yll_spread),[100 T0+1]))
    
    yld_spread = squeeze(sum(model_run.yld_1yr,1)); % years living with the disease
    % Wikipedia: YLD = I x DW x L, where I = number of incident cases in the population, DW = disability weight of specific condition, and L = average duration of the case until remission or death (years)
    % I x L = P therefore YLD = P x DW, which is done in the model
    % model_run.yld_1yr is a 2 x 100 x num_year_1980_2100 matrix 
    % sum over genders to get a 100 x 101 matrix of 100 age groups (0 to 99) versus 101 years (2000 to 2100)
    assert(isequal(size(yld_spread),[100 T0+1]))

    outmat = yll_spread + yld_spread;
    % DALY = YLL + YLD
end % end function make_daly_mat





function out = make_coverage_vec(start_year,num_year_divisions,dt,end_year,yleft,xright,yright)
% yleft contains Montagu coverage values from 1980 to last_available_year
% xright contains boundary years from (last_available_year+dt) onwards
% yright contains coverage values from (last_available_year+dt) to 2101 for boundary years in xright

% expand yleft
% expand yright
% join them
    assert(yleft(1)==0)
    x_vec_before_1980 = start_year:dt:(1980-dt);
    num_time_steps_before_1980 = length(x_vec_before_1980);
    y_vec_before_1980 = repmat(yleft(1),1,num_time_steps_before_1980);

    x_vec_1980_2019 = 1980:dt:(2019-dt);
    num_time_steps_1980_2019 = length(x_vec_1980_2019);
    y_vec_1980_2019 = repmat(yleft(1:(end-1)),num_year_divisions,1);
    assert(isequal(size(y_vec_1980_2019),[num_year_divisions 2018-1979]))
    y_vec_1980_2019 = y_vec_1980_2019(:);
    assert(isequal(size(y_vec_1980_2019),[num_time_steps_1980_2019 1]))

    assert(xright(1)==2019)
    outright = interp1(xright, yright, 2019:dt:end_year, 'linear');

    out = [y_vec_before_1980 y_vec_1980_2019' outright];
    assert(isequal(size(out),size(start_year:dt:end_year)))

end


