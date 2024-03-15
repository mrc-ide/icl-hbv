function output = HBVmodelGAVI_Montagu(cohort_run, seed_prev_string, ...
    num_states, num_sexes, num_treat_blocks, ...
    num_year_divisions, dt, ages, num_age_steps, num_age_bins, start_year, T0, ...
    first_year_of_vacc, last_year_of_vacc, num_vacc_years, vacc_years_pos_vec, ...
    theta, ECofactor, treat_start_year, treat_coverage_in_2016_CDA, ...
    demog, p_ChronicCarriage, Prog, Transactions, treat_cov_in_eligible_2016_factual)


%%% demog = parametes to do with background demography and fitted parts of
%%% the model.
%%%

%% Establish basic simulation parameters
agegroups_10yr = 1 + floor(ages / 10); % categorises the ages into age-groups of 10 year width; 1 x 1000 double; [1 1 ... 10 10], each number present 100 times
agegroups_5yr = 1 + floor(ages / 5); % categorises the ages into age-groups of 5 year width; 1 x 1000 double; [1 1 ... 20 20], each number present 50 times
agegroups_1yr = 1 + floor(ages); % categorises the ages into age-groups of 1 year width; 1 x 1000 double; [1 1 ... 100 100], each number present 10 times
					
i6mo = find(ages >= 0.5, 1); % markers for key age boundaries
i1y = find(ages >= 1, 1);
%i2y = find(ages >= 2, 1);
%i3y = find(ages >= 3, 1);
%i4y = find(ages >= 4, 1);
i5y = find(ages >= 5, 1);
i12y = find(ages >= 12, 1);
i15y = find(ages >= 15, 1);

end_year = start_year + T0; % 2101
TimeSteps = start_year:dt:end_year; % 1 x 2101 double; [1890 1890.1 1890.2 ... 2099.8 2099.9 2100 2100.1 ... 2100.8 2100.9 2101]


%% Establish intervention parameters
% Definition of natural history model (number of stages)
% Definition of stocks (first dimension is age, second dimension is status)


% ---- Intervention Parameters ----

Efficacy_BirthDoseVacc_HbEAg = demog.Efficacy_BirthDoseVacc_HbEAg;
% the efficacy of the vaccination i.e. the proportion of all treated e+ mothers that do not infect their babies
Efficacy_BirthDoseVacc_HbSAg = demog.Efficacy_BirthDoseVacc_HbSAg;
Efficacy_Treatment = 0.98;
Efficacy_InfantVacc = demog.Efficacy_InfantVacc;

p_VerticalTransmission_HbSAg_NoIntv = demog.p_VerticalTransmission_HbSAg_NoIntv; 
% probability of transmission from an HBeAg-, HBsAg+ mother to her baby without intervention
p_VerticalTransmission_HbSAg_BirthDoseVacc = p_VerticalTransmission_HbSAg_NoIntv * (1 - Efficacy_BirthDoseVacc_HbSAg);

p_VerticalTransmission_HbEAg_NoIntv = demog.p_VerticalTransmission_HbEAg_NoIntv;
p_VerticalTransmission_HbEAg_BirthDoseVacc = p_VerticalTransmission_HbEAg_NoIntv * (1 - Efficacy_BirthDoseVacc_HbEAg); 
% probability of transmission from an HBeAg+ mother to her baby after the baby is given BD vaccination

p_VerticalTransmission_Tr_NoIntv = p_VerticalTransmission_HbEAg_NoIntv * (1 - Efficacy_Treatment); 
% probability of transmission from an HBeAg+ mother on treatment to her baby without intervention
p_VerticalTransmission_Tr_BirthDoseVacc = 0.005;



% ---- Load Epidemiological Parameters from fitting procedure ----
beta_U5 = demog.beta_U5;                                      % Rate of horizontal transmission between susceptible and infected persons - UNDER FIVE
beta_1to15 = demog.beta_1to15;                                % Rate of generation transmission between susceptible and infected persons - All Ages
beta_5plus = demog.beta_5plus;
ReducInTransmission = demog.ReducInTransmission;              % Fractional reduction in transmission
YearReducInTransmission = demog.YearReducInTransmission;      % Turning point year for reduction
DurReducInTransmission = 15;                                  % Time taken to complete change.
PriorTDFTreatRate = demog.PriorTDFTreatRate;                  % Treatment rate since 2005
PriorLAMTreatRate = demog.PriorLAMTreatRate;                  % Treatment rate since 2005



% ----- Infection-related parameters -----
beta_scaler = ReducInTransmission ./ (1 + exp( (TimeSteps - (YearReducInTransmission)) ./ (DurReducInTransmission / 10) ));

beta_U5_SAg_itt = beta_U5 * (1 - ReducInTransmission) + beta_U5 * zeros(size(beta_scaler));  % NOT TIME DEPENDENT
beta_U5_EAg_itt = min(1.0, beta_U5_SAg_itt * ECofactor);

beta_1to15_SAg_itt = beta_1to15 * (1 - ReducInTransmission) + beta_1to15 * beta_scaler;  % TIME DEPENDENT
beta_1to15_EAg_itt = min(1.0, beta_1to15_SAg_itt * ECofactor);

beta_5plus_SAg_itt = beta_5plus * (1 - ReducInTransmission) + beta_5plus * beta_scaler;  % TIME DEPENDENT
beta_5plus_EAg_itt = min(1.0, beta_5plus_SAg_itt * ECofactor);



% pop_mat-stocks are (infec, age, sex(1 = women, 2 = men), accessible*)   {*accessible
% specifies whether this person can be reached by treatment progs, 1 = no, 2 = yes}



% Demography
StartPop = demog.Pop_byAgeGroups_1950(agegroups_1yr, :) * dt;
% demog.Pop_byAgeGroups_1950 is a 100 x 2 double of 100 age groups and 2 sexes
% agegroups_1yr is a 1 x 1000 double; [1 1 ... 100 100], each number present 10 times
% start population size like 1950 population
% expanding demog.Pop_byAgeGroups_1950 from 1 year age steps to 0.1 year age steps
% each age group repeated 10 times therefore divide each entry by 10

pop_mat = zeros(num_states, num_age_steps, num_sexes, num_treat_blocks);
% dimensions: disease states, age, sex, treatment

if strcmp(seed_prev_string, 'Cui')
    StartPrev_byAgeGroups = demog.HBsAg_prevs_middle_year_1;
    assert(isequal(size(StartPrev_byAgeGroups), [18 num_sexes]))
    if any(isnan(StartPrev_byAgeGroups))
       nan_positions = isnan(StartPrev_byAgeGroups);
       first_non_nan_pos = find(~nan_positions(:, 1), 1);
       StartPrev_byAgeGroups(1:first_non_nan_pos, :) = repmat(StartPrev_byAgeGroups(first_non_nan_pos, :), first_non_nan_pos, 1);
       last_non_nan_pos = find(~nan_positions(:, 1), 1, 'last'); % find the last non-NaN position
       StartPrev_byAgeGroups(last_non_nan_pos:end, :) = repmat(StartPrev_byAgeGroups(last_non_nan_pos, :), size(nan_positions, 1) - last_non_nan_pos + 1, 1);
    end
    StartPrev_byAgeGroups = [ ...
        StartPrev_byAgeGroups(1:(end - 1), :); ...
        repmat(StartPrev_byAgeGroups(end, :), 3, 1) ...
        ];
    % demog.HBsAg_prevs_middle_year_1 is a 18 x 2 double of age group by sex
    % demog.HBsAg_prevs_middle_year_1 age groups: 0--4 5--9 10--14 15--19 20--24 25--29 30--34 35--39 40--44 45--49 50--54 55--59 60--64 65--69 70--74 75--79 80--84 85+
    assert(isequal(size(StartPrev_byAgeGroups), [20 num_sexes]))

    NumSAg = StartPrev_byAgeGroups(agegroups_5yr, :) .* StartPop;
    % agegroups_5yr is a 1 x 1000 double; [1 1 ... 20 20], each number present 50 times
    % expanding StartPrev_byAgeGroups from 5 year age steps to 0.1 year age steps
    NumNotSAg = (1 - StartPrev_byAgeGroups(agegroups_5yr, :)) .* StartPop;
elseif strcmp(seed_prev_string, 'CDA')
    under_6_pos_vec_len = length(find(ages < 6.0));
    over_6_pos_vec_len = length(find(ages >= 6.0));
    StartPrev_byAgeGroups = [ ...
        repmat(demog.country_HBsAg_prevalences_by_ages_mid_1_young_old(1), under_6_pos_vec_len, num_sexes); ...
        repmat(demog.country_HBsAg_prevalences_by_ages_mid_1_young_old(2), over_6_pos_vec_len, num_sexes) ...
        ];
    % apply prevalence in 5-year-olds to 0 to 6 year olds; apply prevalence in all ages to 6 to 99 year olds
    assert(isequal(size(StartPrev_byAgeGroups), [num_age_steps num_sexes]))

    NumSAg = StartPrev_byAgeGroups .* StartPop;
    NumNotSAg = (1 - StartPrev_byAgeGroups) .* StartPop;
elseif strcmp(seed_prev_string, 'WHO')
    under_5_pos_vec_len = length(find(ages < 5.0)); % ages 0.0 to 4.9 together occupy the first 50 spaces in the ages vector
    over_5_pos_vec_len = length(find(ages >= 5.0)); % ages 5.0 to 99.9 together occupy the last 950 spaces in the ages vector
    assert(under_5_pos_vec_len + over_5_pos_vec_len == num_age_steps)
    StartPrev_byAgeGroups = [ ...
        repmat(demog.country_HBsAg_prevalences_by_ages_prevacc_young_old(1), under_5_pos_vec_len, num_sexes); ...
        repmat(demog.country_HBsAg_prevalences_by_ages_prevacc_young_old(2), over_5_pos_vec_len, num_sexes) ...
        ];
    assert(isequal(size(StartPrev_byAgeGroups), [num_age_steps num_sexes]))

    NumSAg = StartPrev_byAgeGroups .* StartPop;
    NumNotSAg = (1 - StartPrev_byAgeGroups) .* StartPop;
else
    assert(strcmp(seed_prev_string, 'Schweitzer'))

    StartPrev_byAgeGroups = repmat(demog.country_HBsAg_prevalences_by_ages_mid_1_young_old, num_age_steps, num_sexes);
    assert(isequal(size(StartPrev_byAgeGroups), [num_age_steps num_sexes]))

    NumSAg = StartPrev_byAgeGroups .* StartPop;
    NumNotSAg = (1 - StartPrev_byAgeGroups) .* StartPop;
end    

pop_mat(1, :, :, 1) = NumNotSAg;
pop_mat(3, :, :, 1) = 0.5 * NumSAg;
pop_mat(4, :, :, 1) = 0.5 * NumSAg;

 
cov_BirthDose_itt = expand_vacc_cov_vec(demog.BirthDose(:, vacc_years_pos_vec), first_year_of_vacc, last_year_of_vacc, num_vacc_years, num_year_divisions, dt, start_year, end_year, cohort_run);
cov_InfantVacc_0_itt = expand_vacc_cov_vec(demog.InfantVacc_0(:, vacc_years_pos_vec), first_year_of_vacc, last_year_of_vacc, num_vacc_years, num_year_divisions, dt, start_year, end_year, cohort_run);
%cov_InfantVacc_1_itt = expand_vacc_cov_vec(demog.InfantVacc_1, first_year_of_vacc, last_year_of_vacc, num_vacc_years, num_year_divisions, dt, start_year, end_year, cohort_run);
%cov_InfantVacc_2_itt = expand_vacc_cov_vec(demog.InfantVacc_2, first_year_of_vacc, last_year_of_vacc, num_vacc_years, num_year_divisions, dt, start_year, end_year, cohort_run);
%cov_InfantVacc_3_itt = expand_vacc_cov_vec(demog.InfantVacc_3, first_year_of_vacc, last_year_of_vacc, num_vacc_years, num_year_divisions, dt, start_year, end_year, cohort_run);
assert(isequal(size(cov_BirthDose_itt), size(TimeSteps)))
assert(isequal(size(cov_InfantVacc_0_itt), size(TimeSteps)))
%assert(isequal(size(cov_InfantVacc_1_itt), size(TimeSteps)))
%assert(isequal(size(cov_InfantVacc_2_itt), size(TimeSteps)))
%assert(isequal(size(cov_InfantVacc_3_itt), size(TimeSteps)))



assert(isequal(size(Prog), size(zeros(num_states, num_states)))); % Non-Age Specific Prog parameters stored as (from, to)




% Prepare for simulation

% prepare storage containers, for outputs once per year
% breakdowns by age/sex
[Tot_Pop_5yr, Deaths_Severe_Acute_5yr, Deaths_Comp_Cirr_5yr, Deaths_Decomp_Cirr_5yr, Deaths_Liver_Cancer_5yr, PrevSAg_5yr, NumSAg_5yr] ...
    = deal(-99 * ones(num_sexes, max(agegroups_5yr), T0 + 1)); % 1890 + T0 gives T0 + 1 years in total
[Tot_Pop_1yr, Tot_Pop_and_dead_1yr, Prev_Susceptible_1yr, Incid_blue_circle_1yr, Incid_infections_all_1yr, Incid_Severe_Acute_1yr, Incid_Non_Severe_Acute_1yr, Prev_Non_Severe_Acute_1yr, Incid_chronic_horizontal_1yr, Incid_chronic_all_1yr, Prev_Immune_Tolerant_1yr, Incid_Immune_Reactive_1yr, Prev_Immune_Reactive_1yr, Prev_Asymptomatic_Carrier_1yr, Prev_Chronic_Hep_B_1yr, Prev_Comp_Cirr_1yr, Incid_Decomp_Cirr_1yr, Prev_Decomp_Cirr_1yr, Incid_Liver_Cancer_1yr, Incid_Catchup_1yr, Incid_Recovered_State4_1yr, Incid_Recovered_State14_1yr, Incid_Recovered_State15_1yr, Incid_TDF_treat_1yr, Prev_TDF_treat_1yr, Prev_3TC_Treatment_1yr, Prev_Failed_3TC_treat_1yr, Deaths_Comp_Cirr_1yr, Deaths_Decomp_Cirr_1yr, Deaths_Liver_Cancer_1yr, Deaths_TDF_treat_1yr, Deaths_Failed_treat_1yr, Deaths_Severe_Acute_1yr, NumSAg_1yr, PrevSAg_1yr, NumEAg_chronic_acute_1yr, yld_1yr, Incid_Deaths_1yr, Prev_Deaths_1yr] ...
    = deal(-99 * ones(num_sexes, max(agegroups_1yr), T0 + 1));
[...
    Incid_infections_all_1yr_approx, Incid_Severe_Acute_1yr_approx, ...
    Incid_Immune_Reactive_1yr_approx, Incid_Decomp_Cirr_1yr_approx, Incid_TDF_treat_1yr_approx, ...
    Incid_Deaths_1yr_approx...
    ] = deal(-99 * ones(num_sexes, max(agegroups_1yr), T0 + 1));
[NumSAg_10yr] = deal(-99 * ones(num_sexes, max(agegroups_10yr), T0 + 1));
							  
% single output per year
[Time, num_births_1yr, num_at_risk_births_1yr, num_births_Chronic_1yr, num_births_NotChronic_1yr, Incid_babies_chronic_1yr, num_eligible_for_infant_vacc_1yr, num_six_month_olds_1yr, proportion_HBeAg_births_1yr] = deal(-99 * ones(1, T0 + 1));

% double output per year

% ----- Simulation -----

itt = 1; % itt increase every time i.e. every 0.1 years; goes from 1 to 2101 (length of TimeSteps)
OutputEventNum = 1; % OutputEventNum increase every year; goes from 1 to 212
moving_to_treatment = zeros(size(pop_mat));
ReCalibAt1950Done = false;
initiated_treatment = false;
CapturePopMat2000Done = false;
CapturePopMat2005Done = false;
CaptureFOI2000Done = false;
[...
    NewInfections, SevereAcute, NonsevereAcute, moving_btw_states, ...
    into_blue_node, into_state_8, ...
    into_state_10, ...
    state_6_to_11, state_7_to_11, state_8_to_11 ...
    ] = deal(zeros(1, num_age_steps, num_sexes, num_treat_blocks));
into_babies_infected = zeros(num_sexes, 1);
into_babies_chronic = zeros(num_sexes, 1);
babies_infected = 0;
babies_ChronicCarriage = 0;
num_babies_Infected = 0;
num_babies_Chronic = 0;
num_babies_NotChronic = 0;
female_multiplier = 0;
male_multiplier = 0;


for time = TimeSteps 
       
    
    % Update mortality and fertility rates
    mu = zeros(num_states, num_age_steps, num_sexes, num_treat_blocks);
    if cohort_run && time > 2101 % if cohort_run == true and time > 2101
        assert(isequal(size(demog.MortalityRate_Women), size(demog.MortalityRate_Men)))
        assert(isequal(size(demog.MortalityRate_Women), [212 num_age_steps]))
        % demog.MortalityRate_Women is a 212 x 21 matrix of mortality rates of 21 age groups (0--0, 1--4, 5--9, 10--14, ..., 90--94, 95--99) for every year from 1890 to 2101
        % OutputEventNum ranges from 1 to 212
        % agegroups_1yr = [1 1 1 ... 100 100 100], each number 10 times
        % selected vector copied across disease states and treatments
        % use rates for 2101 if time > 2101
        mu(:, :, 1, :) = repmat(demog.MortalityRate_Women(end, :), [num_states 1 num_treat_blocks]);
        mu(:, :, 2, :) = repmat(demog.MortalityRate_Men(end, :), [num_states 1 num_treat_blocks]);
        % demog.MortalityRate_Women is a 212 x 1000 matrix of mortality rates of 1000 age groups (ages 0, 0.1, 0.2,..., 99.8, 99.9) for every year from 1890 to 2101
        assert(isequal(size(demog.fert), [num_age_steps 212]))
        fert = demog.fert(:, end);
        % demog.fert is a 1000 x 212 matrix; ages in 0.1 year jumps versus 212 years
        assert(isequal(size(fert), [num_age_steps 1]))
        assert(length(demog.net_migration) == 212)
        net_migration = demog.net_migration(end);
        assert(length(demog.sex_ratios) == 212)
        sex_ratio = demog.sex_ratios(end);
    else % if cohort_run == false i.e. doing a calendar year run, or if cohort_run == true and time <= 2101
        assert(~cohort_run | time <= 2101)
        % demog.MortalityRate_Women is a 212 x 21 matrix of mortality rates of 21 age groups (0--0, 1--4, 5--9, 10--14, ..., 90--94, 95--99) for every year from 1890 to 2101
        % OutputEventNum ranges from 1 to 212
        % agegroups_1yr = [1 1 1 ... 100 100 100], each number 10 times
        % selected vector copied across disease states and treatments
        mu(:, :, 1, :) = repmat(demog.MortalityRate_Women(OutputEventNum, :), [num_states 1 num_treat_blocks]);
        mu(:, :, 2, :) = repmat(demog.MortalityRate_Men(OutputEventNum, :), [num_states 1 num_treat_blocks]);
        % demog.MortalityRate_Women is a 212 x 1000 matrix of mortality rates of 1000 age groups (ages 0, 0.1, 0.2,..., 99.8, 99.9) for every year from 1890 to 2101
        fert = demog.fert(:, OutputEventNum);
        % demog.fert is a 1000 x 212 matrix; ages in 0.1 year jumps versus 212 years
        assert(isequal(size(fert), [num_age_steps 1]))
        assert(length(demog.net_migration) == 212)
        net_migration = demog.net_migration(OutputEventNum);
        sex_ratio = demog.sex_ratios(OutputEventNum);
    end
    assert(length(net_migration) == 1)
    net_migration = repmat(net_migration, [num_states num_age_steps num_sexes num_treat_blocks]);


    
    % Compute Outputs once per year
    if rem(time, 1) == 0 
    % remainder after division by 1 is zero i.e. one has a whole number
    % only saves variables in this "for" loop every 10 time steps (or once a year, since dt=0.1)
    % contains variables calculated from variables that have been accumulated over the previous year; hence "OutputEventNum - 1"


        Time(OutputEventNum) = time;
        % OutputEventNum ranges from 1 to 212

    
        for k = 1:num_sexes

                
            for ag = 1:max(agegroups_5yr) % 1:20
														
                if OutputEventNum > 1
                    % incidences and deaths are accummulated every 0.1 years in the Disease Progression loop and 
                    % then assigned to a particular year at the beginning of that year, after which they are zeroed
                    
                    Deaths_Comp_Cirr_5yr(k, ag, OutputEventNum - 1) = squeeze(sum(sum(state_6_to_11(1, agegroups_5yr == ag, k, :), 2), 4));

                    Deaths_Decomp_Cirr_5yr(k, ag, OutputEventNum - 1) = squeeze(sum(sum(state_7_to_11(1, agegroups_5yr == ag, k, :), 2), 4));

                    Deaths_Liver_Cancer_5yr(k, ag, OutputEventNum - 1) = squeeze(sum(sum(state_8_to_11(1, agegroups_5yr == ag, k, :), 2), 4));

                end
                
            end % end agegroups_5yr for loop

            
            for ag = 1:max(agegroups_1yr) % 1:100
 
                if OutputEventNum > 1

                    Incid_Liver_Cancer_1yr(k, ag, OutputEventNum - 1) = squeeze(sum(sum(into_state_8(1, agegroups_1yr == ag, k, :), 2), 4));

                end
                
            end % end agegroups_1yr for loop
                
             
        end % end sexes for loop

        
        [into_blue_node, into_state_8, ...
            into_state_10, ...
            state_6_to_11, state_7_to_11, state_8_to_11...
            ] = deal(zeros(1, num_age_steps, num_sexes, num_treat_blocks));
        into_babies_infected = zeros(num_sexes, 1);    
        into_babies_chronic = zeros(num_sexes, 1);    
        num_babies_Infected = 0;
        num_babies_Chronic = 0;
        num_babies_NotChronic = 0;
        
    end % end "rem(time, 1) == 0" if statement
    
    
    % Horizontal Transmission (infection and Chronic Carriage)
    
    FOI = zeros(1, num_age_steps, num_sexes, num_treat_blocks);
    
    % i: Transmission Between Children (1y--5y olds)
    FOI(1, i1y:(i5y - 1), :, :) = ...
        beta_U5_SAg_itt(itt) * sum(sum(sum(sum(pop_mat([4:8 13], i1y:(i5y - 1), :, :))))) / sum(sum(sum(sum(pop_mat([1:10 12:15], i1y:(i5y - 1), :, :))))) ...
        + beta_U5_EAg_itt(itt) * sum(sum(sum(sum(pop_mat([2:3 14:15], i1y:(i5y - 1), :, :))))) / sum(sum(sum(sum(pop_mat([1:10 12:15], i1y:(i5y - 1), :, :)))));
    
    % ii. Transmission between 1--15 year olds
    FOI(1, i1y:(i15y - 1), :, :) = FOI(1, i1y:(i15y - 1), :, :) + ...
        beta_1to15_SAg_itt(itt) * sum(sum(sum(sum(pop_mat([4:8 13], i1y:(i15y - 1), :, :))))) / sum(sum(sum(sum(pop_mat([1:10 12:15], i1y:(i15y - 1), :, :))))) ...
        + beta_1to15_EAg_itt(itt) * sum(sum(sum(sum(pop_mat([2:3 14:15], i1y:(i15y - 1), :, :))))) / sum(sum(sum(sum(pop_mat([1:10 12:15], i1y:(i15y - 1), :, :)))));
    
    % iii.: Transmission Between 5+ and Adults (Assuming equal risks
    % for all persons 5y--100y)
    FOI(1, i5y:end, :, :) = FOI(1, i5y:end, :, :) + ...
        beta_5plus_SAg_itt(itt) * sum(sum(sum(sum(pop_mat([4:8 13], i5y:end, :, :))))) / sum(sum(sum(sum(pop_mat([1:10 12:15], i5y:end, :, :))))) ...
        + beta_5plus_EAg_itt(itt) * sum(sum(sum(sum(pop_mat([2:3 14:15], i5y:end, :, :))))) / sum(sum(sum(sum(pop_mat([1:10 12:15], i5y:end, :, :)))));
    
    
    % Disease Progression
    next_pop_mat = pop_mat;
    for tr = 1:length(Transactions.From)
        transaction_vals = Transactions.Values{tr};
        transaction_vals = transaction_vals(:);
        %assert(all(transaction_vals >= 0))
        %assert(all(transaction_vals <= 1))
        moving_btw_states(1, :, :, :) = pop_mat(Transactions.From(tr), :, :, :) .* Transactions.Values{tr};
        next_pop_mat(Transactions.From(tr), :, :, :) = next_pop_mat(Transactions.From(tr), :, :, :) + dt * ( -moving_btw_states ); % move people out of "from" state
        next_pop_mat(Transactions.To(tr), :, :, :) = next_pop_mat(Transactions.To(tr), :, :, :)  + dt * ( +moving_btw_states ); % move people into "to" state

        if (Transactions.To(tr) == 8)
            into_state_8 = into_state_8 + dt * moving_btw_states;
        end

        if ((Transactions.From(tr) == 6) && (Transactions.To(tr) == 11))
            state_6_to_11 = state_6_to_11 + dt * moving_btw_states;
        end

        if ((Transactions.From(tr) == 7) && (Transactions.To(tr) == 11))
            state_7_to_11 = state_7_to_11 + dt * moving_btw_states;
        end

        if ((Transactions.From(tr) == 8) && (Transactions.To(tr) == 11))
            state_8_to_11 = state_8_to_11 + dt * moving_btw_states;
        end

    end % end Disease Progression for loop
    % multiply by dt because quantity is calculated every 0.1 years therefore needs to be divided by 10
    


    
    assert(squeeze(sum(sum(sum(sum(pop_mat([12 13], :, :, 1), 1), 2), 3), 4)) == 0)
    % (Time-dependent) Baseline Transition to TDF-Treatment
    if (time >= treat_start_year)
    % 2016 must be the first year with nonzero treatment 

        if ~initiated_treatment
            if ~exist('treat_cov_in_eligible_2016_factual', 'var')
                num_in_treatment = sum(sum(sum(sum(pop_mat(10, :, :, :), 1), 2), 3), 4);
                assert(num_in_treatment == 0) % no one is in treatment
                prev_pop = sum(sum(sum(sum(pop_mat([2:8 10 12:15], :, :, :), 1), 2), 3), 4); 

                total_num_to_move_to_treat = treat_coverage_in_2016_CDA * prev_pop;
                eligible_pop = sum(sum(sum(sum(pop_mat([3 5:7 10], i12y:end, :, :), 1), 2), 3), 4); % only those over 12 years old allowed into treatment
                assert(total_num_to_move_to_treat < eligible_pop)
                treat_cov_in_eligible_2016_factual = total_num_to_move_to_treat / eligible_pop;
            end
            next_pop_mat([3 5 6 7], i12y:end, :, :) = next_pop_mat([3 5 6 7], i12y:end, :, :) - pop_mat([3 5 6 7], i12y:end, :, :) * treat_cov_in_eligible_2016_factual;
            % Every compartment in the eligible-for-treatment states in next_pop_mat must have a number subtracted from it 
            % such that the total number subtracted from the eligible-for-treatment states is in_treatment_2016
            % i.e. in_treatment_2016 = sum(sum(sum(sum(pop_mat([3 5 6 7], :, :, :), 1), 2), 3), 4) * treat_cov_in_eligible_2016_factual = sum(sum(sum(sum(pop_mat([3 5 6 7], :, :, :) * treat_cov_in_eligible_2016_factual, 1), 2), 3), 4)
            % Hence, treat_cov_in_eligible_2016_factual scales each compartment in pop_mat([3 5 6 7], :, :, :) such that pop_mat([3 5 6 7], :, :, :) * treat_cov_in_eligible_2016_factual subtracts the same proportion of people from each compartment in each of the eligible-for-treatment states in order to subtract a total of in_treatment_2016 from the eligible-for-treatment states.
            next_pop_mat(10, i12y:end, :, :) = next_pop_mat(10, i12y:end, :, :) + sum(pop_mat([3 5 6 7], i12y:end, :, :) * treat_cov_in_eligible_2016_factual, 1);

            assert(sum(sum(sum(sum(next_pop_mat(10, 1:(i12y - 1), :, :), 1), 2), 3), 4) == 0) % ensure that no 0 to 11 year olds are in treatment
            num_in_treatment = sum(sum(sum(sum(next_pop_mat(10, i12y:end, :, :), 1), 2), 3), 4);
            eligible_pop = sum(sum(sum(sum(pop_mat([3 5:7 10], i12y:end, :, :), 1), 2), 3), 4); 
            assert(num_in_treatment / eligible_pop >= treat_coverage_in_2016_CDA)
            % treatment coverage amongst treatment-eligible people will be greater than treatment coverage amongst HBsAg+ people, except if treatment coverage is 0
            treat_coverage_2016 = num_in_treatment / eligible_pop;
            assert(abs(treat_coverage_2016 - treat_cov_in_eligible_2016_factual) < 1e-10)

            initiated_treatment = true;
        else
            assert(initiated_treatment) % ensure that, each time this code is encountered, treatment has already been initiated

            assert(demog.PriorTDFTreatRate >= 0)
            moving_to_treatment([3 5 6 7], i12y:end, :, :) = pop_mat([3 5 6 7], i12y:end, :, :) .* demog.PriorTDFTreatRate;
            next_pop_mat([3 5 6 7], i12y:end, :, :) = next_pop_mat([3 5 6 7], i12y:end, :, :) + dt * ( -moving_to_treatment([3 5 6 7], i12y:end, :, :) );
            next_pop_mat(10, i12y:end, :, :) = next_pop_mat(10, i12y:end, :, :) + dt * ( +sum(moving_to_treatment(:, i12y:end, :, :), 1) );
            assert(max(moving_to_treatment(:)) >= 0)
            into_state_10 = into_state_10 + dt * ( +sum(moving_to_treatment, 1) );
            assert(sum(sum(sum(sum(moving_to_treatment(:, 1:(i12y - 1), :, :), 1), 2), 3), 4) == 0) % ensure that no 0 to 11 year olds are moved to treatment
            if demog.PriorTDFTreatRate > 0
                assert(any(into_state_10(:) > 0))
            end
        end
    end % end treatment if statement
    

    
    
    % Infection process
    NewInfections = pop_mat(1, :, :, :) .* FOI;
    % NewInfections is the blue node
    % number of susceptibles times FOI i.e. number of new infections within population, excluding babies (since FOI is 0 for babies)
    % 1 x num_age_steps x 2 x 2 double i.e. 1 x 1000 x 2 x 2 double
    %NewChronicCarriage = NewInfections .* p_ChronicCarriage;
    SevereAcute = NewInfections * theta;
    NonsevereAcute = NewInfections - SevereAcute;
    
    % Transitions dependent on the blue node, which does not have a number and is therefore not in Prog or Transactions
    % multiply by dt, since FOI is an annual rate
    into_blue_node = into_blue_node + dt * NewInfections;
    next_pop_mat(1, :, :, :) = next_pop_mat(1, :, :, :) + dt * ( -NewInfections );
    next_pop_mat(14, :, :, :) = next_pop_mat(14, :, :, :) + dt * ( +NonsevereAcute );
    next_pop_mat(15, :, :, :) = next_pop_mat(15, :, :, :) + dt * ( +SevereAcute );
    
    
    % Infanct vaccination occurring at exactly six months
    % do not multiply by dt, since one is vaccinating cov_InfantVacc_itt(itt)% of people in next_pop_mat(1, i6mo, :, :), 
    % after which this cohort ages and moves to the next age bin
    % if divides all babies born in a year into 10 groups and vaccinatates cov_InfantVacc_itt(itt)% of each group, 
    % then one will have vaccinated cov_InfantVacc_itt(itt)% of all babies born in that year
    transfer_to_vacc = cov_InfantVacc_0_itt(itt) * next_pop_mat(1, i6mo, :, :) * Efficacy_InfantVacc; % the 0.95 represent a take-type vaccine efficacy of 95%.
    next_pop_mat(1, i6mo, :, :) = next_pop_mat(1, i6mo, :, :) - transfer_to_vacc;
    next_pop_mat(9, i6mo, :, :) = next_pop_mat(9, i6mo, :, :) + transfer_to_vacc;
    
        
    
    % Natural Mortality
    % Do not apply background mortality to the HBV deaths state, since we want people in all countries to be treated as if they would have lived until 84 if they had not died of HBV. This is done outside of the model in the main script. 
    mu(11, :, :, :) = 0.0;
    next_pop_mat = next_pop_mat + dt * ( -next_pop_mat .* mu );
    % in the "for time=TimeSteps", which runs 10 times per year, therefore divide effect of mu by 10
    net_migration(11, :, :, :) = 0.0;
    next_pop_mat = next_pop_mat + dt * ( +next_pop_mat .* net_migration );
    % "+" because net_migration = (number of immigrants - number of emigrants) / population size
    
    % Update Stocks
    pop_mat = next_pop_mat;
    clear next_pop_mat
    
    % now age everyone
    pop_mat(:, 2:num_age_steps, :, :) = pop_mat(:, 1:(num_age_steps - 1), :, :);
    pop_mat(:, 1, :, :) = 0; % set number of babies to 0 (babies will be born next)
    
    
    % fill-out with new births in this time-step:
    births_toNonInfectiousWomen = sum( fert' .* sum(sum(pop_mat([1 9], :, 1, :), 1), 4) ); % Susecptible, Immune
    births_toHbEAgWomen = sum(fert' .* sum(sum(pop_mat([2:3 14:15], :, 1, :), 1), 4)); % Immune Tolerant, Immune Reactive
    births_toHbSAgWomen = sum(fert' .* sum(sum(pop_mat([4:8 13], :, 1, :), 1), 4)); % All other stages (other infected women)
    births_toTrWomen = sum(fert' .* sum(sum(pop_mat([10 12], :, 1, :), 1), 4)); % Women on Treatment
    
    
    births_toInfectiousWomen = births_toHbEAgWomen + births_toHbSAgWomen + births_toTrWomen;
    births_Total = births_toNonInfectiousWomen + births_toInfectiousWomen;
    assert(length(births_Total) == 1)
    

    babies_infected = ... % a 1 x 1 double
        ...
        births_toHbSAgWomen * (1 - cov_BirthDose_itt(itt)) * p_VerticalTransmission_HbSAg_NoIntv ...
        + births_toHbSAgWomen * cov_BirthDose_itt(itt) * p_VerticalTransmission_HbSAg_BirthDoseVacc ...
        ...
        + births_toHbEAgWomen * (1 - cov_BirthDose_itt(itt)) * p_VerticalTransmission_HbEAg_NoIntv ...
        + births_toHbEAgWomen * cov_BirthDose_itt(itt) * p_VerticalTransmission_HbEAg_BirthDoseVacc ...
        ...
        + births_toTrWomen * (1 - cov_BirthDose_itt(itt)) * p_VerticalTransmission_Tr_NoIntv ...
        + births_toTrWomen * cov_BirthDose_itt(itt) * p_VerticalTransmission_Tr_BirthDoseVacc;
    babies_ChronicCarriage = p_ChronicCarriage(1, 1, 1, 1) * babies_infected;
    
    babies_NotChronicCarriage = births_Total - babies_ChronicCarriage;
    assert(length(babies_infected) == 1)
    assert(length(babies_ChronicCarriage) == 1)
    assert(length(babies_NotChronicCarriage) == 1)
    num_babies_Infected = num_babies_Infected + dt * babies_infected;
    num_babies_Chronic = num_babies_Chronic + dt * babies_ChronicCarriage;
    num_babies_NotChronic = num_babies_NotChronic + dt * babies_NotChronicCarriage;

    
    female_multiplier = 1 / (1 + sex_ratio);
    male_multiplier = sex_ratio / (1 + sex_ratio);
    % sex_ratio is number of male births per one female birth
    % 1 / (1 + sex_ratio) + sex_ratio / (1 + sex_ratio) = 1, hence total number of babies not changed
    % male births -> 0 => sex_ratio -> 0 => male_multiplier -> 0 / 1 = 0
    % male births -> infinity => sex_ratio -> infinity => male_multiplier -> 1
    % female births -> 0 => sex_ratio -> infinity => female_multiplier -> 0
    % female births -> infinity => sex_ratio -> 0 => female_multiplier -> 1
    
    pop_mat(1, 1, 1, 1) = female_multiplier * dt * babies_NotChronicCarriage;  % Suscpetible babies
    pop_mat(2, 1, 1, 1) = female_multiplier * dt * babies_ChronicCarriage;     % Babies with chronic carriage
    
    pop_mat(1, 1, 2, 1) = male_multiplier * dt * babies_NotChronicCarriage;  % Suscpetible babies
    pop_mat(2, 1, 2, 1) = male_multiplier * dt * babies_ChronicCarriage;     % Babies with chronic carriage
  
    assert(isequal(size(into_babies_chronic), [num_sexes 1]))
    assert(isequal(size(squeeze(pop_mat(2, 1, :, 1))), [num_sexes 1]))
    into_babies_chronic = into_babies_chronic + squeeze(pop_mat(2, 1, :, 1));
    assert(isequal(size(into_babies_infected), [num_sexes 1]))
    into_babies_infected = into_babies_infected + [female_multiplier * dt * babies_infected; male_multiplier * dt * babies_infected];
    
    
    % Compute Outputs once per year
    if rem(time, 1) == 0 
    % remainder after division by 1 is zero i.e. one has a whole number
    % only saves variables in this "for" loop every 10 time steps (or once a year, since dt=0.1)
    % contains variables calculated from variables that capture an instant of time i.e. variables not accumulated over the previous year
 
        
        % Rescale population sizes of each age group and sex
        if (time >= 1950) 
            ModelPop = squeeze(sum(sum(pop_mat([1:10 12:15], :, :, :), 1), 4));
            assert(isequal(size(ModelPop), [num_age_steps num_sexes]))
            % sum over disease state of alive people and treatment; ModelPop is 1000 x 2 i.e. age groups versus sex
            assert(isequal(size(demog.total_pop_female), [101 152]))
            % demog.total_pop_female is a 101 x 152 matrix of 152 years (1950 to 2101) for 101 age groups (0--0, 1--1, 2--2,..., 98--98, 99--99, 100--100)
            assert(isequal(size(demog.total_pop_male), [101 152]))
            if cohort_run && time > 2101 % if cohort_run == true and time > 2101
                col_index = 152;
            else % if cohort_run == false i.e. doing a calendar year run, or if cohort_run == true and time <= 2101
                assert(~cohort_run | time <= 2101)
                col_index = time - 1949;
            end
            MontaguPopFemale = demog.total_pop_female(1:num_age_bins, col_index);
            % only want ages 0--99
            MontaguPopMale = demog.total_pop_male(1:num_age_bins, col_index);
            MontaguPop = [MontaguPopFemale MontaguPopMale];
            assert(isequal(size(MontaguPop), [num_age_bins num_sexes]))
            MontaguPopExpand = MontaguPop(agegroups_1yr, :) * dt;
            % agegroups_1yr is a 1 x 1000 double; [1 1 ... 100 100], each number present 10 times
            % expanding MontaguPop from 1 year age steps to 0.1 year age steps
            % each age group repeated 10 times therefore divide each entry by 10
            assert(isequal(size(MontaguPopExpand), [num_age_steps num_sexes]))
            ScalerMat = MontaguPopExpand ./ ModelPop;
            ScalerMat(isnan(ScalerMat)) = 0;
            ScalerMat(isinf(ScalerMat)) = 0;
            pop_scaler = repmat(reshape(ScalerMat, [1 num_age_steps num_sexes]), [num_states 1 1 num_treat_blocks]);
            % MontaguPopExpand is sizes of the current year's population over 0.1 year age steps; a 1000 x 2 matrix of ages versus sex
            % add an extra dimension and duplicate it for each disease state and treatment method
            pop_mat = pop_mat .* pop_scaler;
            % scale all parts of pop_mat, including dead people i.e. State 11
        end
        
        
        for k = 1:num_sexes
        
        
            for ag = 1:max(agegroups_10yr) % 1:10
														
                NumSAg_10yr(k, ag, OutputEventNum) = squeeze(sum(sum(sum(pop_mat([2:8 10 12:15], agegroups_10yr == ag, k, :), 1), 2), 4));
                
            end % end agegroups_10yr for loop

                
            for ag = 1:max(agegroups_5yr) % 1:20
                
                state_prev_vec = squeeze(sum(sum(pop_mat(:, agegroups_5yr == ag, k, :), 2), 4)); % k is sex
                assert(isequal(size(state_prev_vec), [num_states 1]))
                
                Tot_Pop_5yr(k, ag, OutputEventNum) = sum(state_prev_vec([1:10 12:15]));

                NumSAg_5yr(k, ag, OutputEventNum) = sum(state_prev_vec([2:8 10 12:15]));
                
            end % end agegroups_5yr for loop

            
            for ag = 1:max(agegroups_1yr) % 1:100
                
                state_prev_vec = squeeze(sum(sum(pop_mat(:, agegroups_1yr == ag, k, :), 2), 4)); % k is sex
                assert(isequal(size(state_prev_vec), [num_states 1]))

                Tot_Pop_1yr(k, ag, OutputEventNum) = sum(state_prev_vec([1:10 12:15]));
                
                Prev_Immune_Reactive_1yr(k, ag, OutputEventNum) = state_prev_vec(3);
                
                Prev_Chronic_Hep_B_1yr(k, ag, OutputEventNum) = state_prev_vec(5);

                Prev_Comp_Cirr_1yr(k, ag, OutputEventNum) = state_prev_vec(6);

                Prev_Decomp_Cirr_1yr(k, ag, OutputEventNum) = state_prev_vec(7);
                
                Prev_TDF_treat_1yr(k, ag, OutputEventNum) = state_prev_vec(10);

                NumSAg_1yr(k, ag, OutputEventNum) = sum(state_prev_vec([2:8 10 12:15]));
                
                Incid_Decomp_Cirr_1yr_approx(k, ag, OutputEventNum) = sum(state_prev_vec .* Prog(:, 7));
                
                Incid_Deaths_1yr_approx(k, ag, OutputEventNum) = sum(state_prev_vec .* Prog(:, 11));

                if ((time > 2016) && (ag >= 13)) % in 2016 people are moved en masse into the treatment state, so incidence is only greater than 0 from 2017 onwards
                    Incid_TDF_treat_1yr_approx(k, ag, OutputEventNum) = sum(state_prev_vec([3 5 6 7])) .* demog.PriorTDFTreatRate;
                else
                    Incid_TDF_treat_1yr_approx(k, ag, OutputEventNum) = 0;
                end

                Incid_infections_all_1yr_approx(k, ag, OutputEventNum) = sum(sum(NewInfections(1, agegroups_1yr == ag, k, :), 2), 4);

                Incid_Severe_Acute_1yr_approx(k, ag, OutputEventNum) = sum(sum(SevereAcute(1, agegroups_1yr == ag, k, :), 2), 4);
                
            end % end agegroups_1yr for loop

            if k == 1
                sex_multiplier = female_multiplier;
            else
                sex_multiplier = male_multiplier;
            end

            assert(Incid_infections_all_1yr_approx(k, 1, OutputEventNum) == 0) % 0-year-olds cannot get horizontal chronic infection (see FOI)							
            Incid_infections_all_1yr_approx(k, 1, OutputEventNum) = sex_multiplier * babies_infected;


        end % end sexes for loop

        proportion_HBeAg_births_1yr(OutputEventNum) = births_toHbEAgWomen / births_toInfectiousWomen;

        % Update counter
        OutputEventNum = OutputEventNum + 1;
        % increases every year
        
    end % end "rem(time, 1) == 0" if statement

      
    % increment the pointer
    itt = itt + 1;
    % increases every 0.1 years
    
end % end "time = TimeSteps" for loop

output.Time = Time; % 1 x (T0 + 1)
output.treat_coverage_2016 = treat_coverage_2016;
output.Tot_Pop_5yr = Tot_Pop_5yr; % 2 x 20 x (T0 + 1)
output.Deaths_Comp_Cirr_5yr = Deaths_Comp_Cirr_5yr; % 2 x 20 x (T0 + 1)
output.Deaths_Decomp_Cirr_5yr = Deaths_Decomp_Cirr_5yr; % 2 x 20 x (T0 + 1)
output.Deaths_Liver_Cancer_5yr = Deaths_Liver_Cancer_5yr; % 2 x 20 x (T0 + 1)
output.NumSAg_5yr = NumSAg_5yr; % 2 x 20 x (T0 + 1)
output.Tot_Pop_1yr = Tot_Pop_1yr; % 2 x 100 x (T0 + 1)
output.Incid_infections_all_1yr_approx = Incid_infections_all_1yr_approx; % 2 x 100 x (T0 + 1)
output.Incid_Severe_Acute_1yr_approx = Incid_Severe_Acute_1yr_approx; % 2 x 100 x (T0 + 1)
output.Prev_Immune_Reactive_1yr = Prev_Immune_Reactive_1yr; % 2 x 100 x (T0 + 1)
output.Prev_Chronic_Hep_B_1yr = Prev_Chronic_Hep_B_1yr; % 2 x 100 x (T0 + 1)
output.Prev_Comp_Cirr_1yr = Prev_Comp_Cirr_1yr; % 2 x 100 x (T0 + 1)
output.Incid_Decomp_Cirr_1yr_approx = Incid_Decomp_Cirr_1yr_approx; % 2 x 100 x (T0 + 1)
output.Prev_Decomp_Cirr_1yr = Prev_Decomp_Cirr_1yr; % 2 x 100 x (T0 + 1)
output.Incid_Liver_Cancer_1yr = Incid_Liver_Cancer_1yr; % 2 x 100 x (T0 + 1)
output.Incid_TDF_treat_1yr_approx = Incid_TDF_treat_1yr_approx; % 2 x 100 x (T0 + 1)
output.Prev_TDF_treat_1yr = Prev_TDF_treat_1yr; % 2 x 100 x (T0 + 1)
output.NumSAg_1yr = NumSAg_1yr; % 2 x 100 x (T0 + 1)
output.Incid_Deaths_1yr_approx = Incid_Deaths_1yr_approx; % 2 x 100 x (T0 + 1)
output.proportion_HBeAg_births_1yr = proportion_HBeAg_births_1yr; % 1 x (T0 + 1)

end % end function HBVmodel_PPT




function cov_itt_vec = expand_vacc_cov_vec(cov_vec, first_year_of_vacc, last_year_of_vacc, num_vacc_years, num_year_divisions, dt, start_year_run, end_year_run, cohort_run)

    assert(end_year_run >= first_year_of_vacc) % coverage is from 1980 onwards
    
    assert(isequal(size(cov_vec), [1 num_vacc_years])) % 1980--2100
    if ((end_year_run > last_year_of_vacc) || (cohort_run)) % i.e. 2101 or cohort_run
        cov_vec = cov_vec;
    else
        assert(end_year_run <= last_year_of_vacc)
        index_end_year = end_year_run - first_year_of_vacc + 1;
        cov_vec = cov_vec(:, 1:index_end_year);
    end
    
    zeros_vec = zeros(1, length(start_year_run:dt:(first_year_of_vacc - dt))); % 1890--(1980-dt)
    
    cov_itt_vec = repmat(cov_vec(1:(end - 1)), num_year_divisions, 1); % 1980--2099 -> 1980--2099.9
    if ((end_year_run > last_year_of_vacc) || (cohort_run)) % i.e. 2101 or cohort_run
        assert(isequal(size(cov_itt_vec), [num_year_divisions ((last_year_of_vacc - first_year_of_vacc + 1) - 1)]))
    else
        assert(isequal(size(cov_itt_vec), [num_year_divisions ((end_year_run - first_year_of_vacc + 1) - 1)]))
    end
    cov_itt_vec = cov_itt_vec(:);
    if ((end_year_run > last_year_of_vacc) || (cohort_run)) % i.e. 2101 or cohort_run
        cov_itt_vec = [zeros_vec cov_itt_vec' repmat(cov_vec(end), 1, num_year_divisions * (end_year_run - last_year_of_vacc) + 1)]; %  2100 -> 2100--2101
        % 1890--(1980-dt), 1980--2099.9, 2100--2101
    else
        cov_itt_vec = [zeros_vec cov_itt_vec' cov_vec(end)]; % 1890--(1980-dt), 1980--2099.9, 2100
        % TimeSteps = [1890 1890.1 1890.2 ... 2099.8 2099.9 2100]
    end

    assert(all(cov_itt_vec >= 0))
    assert(all(cov_itt_vec <= 1))

end % end function expand_vacc_cov_vec


