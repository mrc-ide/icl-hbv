function output = HBVmodel(seed_prev_string,...
    num_states,num_year_divisions,dt,ages,num_age_steps,start_year,T0,...
    theta,ECofactor,treat_start_year,treat_coverage_in_2016,demog,p_ChronicCarriage,Prog,Transactions)


%% Establish basic simulation parameters
agegroups_5yr = 1 + floor(ages / 5); % categorises the ages into age-groups of 5 year width; 1 x 1000 double; [1 1 ... 20 20], each number present 50 times
agegroups_1yr = 1 + floor(ages); % categorises the ages into age-groups of 1 year width; 1 x 1000 double; [1 1 ... 100 100], each number present 10 times
					
i6mo = find(ages >= 0.5, 1); % markers for key age boundaries
i1y = find(ages >= 1, 1);
i5y = find(ages >= 5, 1);
i15y = find(ages >= 15, 1);

end_year = start_year + T0;
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

p_VerticalTransmission_HbSAg_NoIntv = demog.p_VerticalTransmission_HbSAg_NoIntv; % probability of transmission from an HBeAg-, HBsAg+ mother to her baby without intervention
p_VerticalTransmission_HbSAg_BirthDoseVacc = p_VerticalTransmission_HbSAg_NoIntv * (1 - Efficacy_BirthDoseVacc_HbSAg);

p_VerticalTransmission_HbEAg_NoIntv = demog.p_VerticalTransmission_HbEAg_NoIntv;
p_VerticalTransmission_HbEAg_BirthDoseVacc = p_VerticalTransmission_HbEAg_NoIntv * (1 - Efficacy_BirthDoseVacc_HbEAg); 
% probability of transmission from an HBeAg+ mother to her baby after the baby is given BD vaccination

p_VerticalTransmission_Tr_NoIntv = p_VerticalTransmission_HbEAg_NoIntv * (1 - Efficacy_Treatment); % probability of transmission from an HBeAg+ mother on treatment to her baby without intervention
p_VerticalTransmission_Tr_BirthDoseVacc = 0.005;



% ---- Load Epidemiological Parameters from fitting procedure ----
beta_U5 = demog.beta_U5; % Rate of horizontal transmission between susceptible and infected persons - UNDER FIVE
beta_1to15 = demog.beta_1to15;                                % Rate of generation transmission between susceptible and infected persons - All Ages
beta_5plus = demog.beta_5plus;
ReducInTransmission = demog.ReducInTransmission;              % Fractional reduction in transmission
YearReducInTransmission = demog.YearReducInTransmission;      % Turning point year for reduction
DurReducInTransmission = 15;                                  % Time taken to complete change
PriorTDFTreatRate = demog.PriorTDFTreatRate;                  % Treatment rate since 2005



% ----- Infection-relate parameters -----
beta_scaler = ReducInTransmission ./ (1 + exp( (TimeSteps - (YearReducInTransmission)) ./ (DurReducInTransmission / 10) ));

beta_U5_SAg_itt = beta_U5 * (1 - ReducInTransmission) + beta_U5 * zeros(size(beta_scaler));  % NOT TIME DEPENDENT
beta_U5_EAg_itt = min(1.0, beta_U5_SAg_itt * ECofactor);

beta_1to15_SAg_itt = beta_1to15 * (1 - ReducInTransmission) + beta_1to15 * beta_scaler;  % TIME DEPENDENT
beta_1to15_EAg_itt = min(1.0, beta_1to15_SAg_itt * ECofactor);

beta_5plus_SAg_itt = beta_5plus * (1 - ReducInTransmission) + beta_5plus * beta_scaler;  % TIME DEPENDENT
beta_5plus_EAg_itt = min(1.0, beta_5plus_SAg_itt * ECofactor);



% X-stocks are (infec, age, sex(1=women, 2=men), accessible*)   {*accessible
% specifies whether this person can be reached by treatment progs, 1=no, 2=yes}



% Demography
StartPop = demog.Pop_byAgeGroups_1950(agegroups_1yr, :) * dt;
% demog.Pop_byAgeGroups_1950 is a 100 x 2 double of 100 age groups and 2 genders
% agegroups_1yr is a 1 x 1000 double; [1 1 ... 100 100], each number present 10 times
% start population size like 1950 population
% expanding demog.Pop_byAgeGroups_1950 from 1 year age steps to 0.1 year age steps
% each age group repeated 10 times therefore divide each entry by 10

X = zeros(num_states, num_age_steps, 2, 2);
% dimensions: disease states, age, gender, accessible to treatment

if strcmp(seed_prev_string,'Cui')
    StartPrev_byAgeGroups = demog.HBsAg_prevs_middle_year_1;
    assert(isequal(size(StartPrev_byAgeGroups),[18 2]))
    if any(isnan(StartPrev_byAgeGroups))
       nan_positions = isnan(StartPrev_byAgeGroups);
       first_non_nan_pos = find(~nan_positions(:,1), 1);
       StartPrev_byAgeGroups(1:first_non_nan_pos,:) = repmat(StartPrev_byAgeGroups(first_non_nan_pos,:),first_non_nan_pos,1);
       last_non_nan_pos = find(~nan_positions(:,1), 1, 'last'); % find the last non NaN position
       StartPrev_byAgeGroups(last_non_nan_pos:end,:) = repmat(StartPrev_byAgeGroups(last_non_nan_pos,:),size(nan_positions,1)-last_non_nan_pos+1,1);
    end
    StartPrev_byAgeGroups = [StartPrev_byAgeGroups(1:end-1,:); repmat(StartPrev_byAgeGroups(end,:),3,1)];
    % demog.HBsAg_prevs_middle_year_1 is a 18 x 2 double of age group by gender
    % demog.HBsAg_prevs_middle_year_1 age groups: 0--4 5--9 10--14 15--19 20--24 25--29 30--34 35--39 40--44 45--49 50--54 55--59 60--64 65--69 70--74 75--79 80--84 85+
    assert(isequal(size(StartPrev_byAgeGroups),[20 2]))

    NumSAg = StartPrev_byAgeGroups(agegroups_5yr, :) .* StartPop;
    % agegroups_5yr is a 1 x 1000 double; [1 1 ... 20 20], each number present 50 times
    % expanding StartPrev_byAgeGroups from 5 year age steps to 0.1 year age steps
    NumNotSAg = (1 - StartPrev_byAgeGroups(agegroups_5yr, :)) .* StartPop;
elseif strcmp(seed_prev_string,'CDA')
    StartPrev_byAgeGroups = [repmat(demog.country_HBsAg_prevalences_by_ages_mid_1_young_old(1),num_year_divisions*(5.9-0.0)+1,2); ...
        repmat(demog.country_HBsAg_prevalences_by_ages_mid_1_young_old(2),num_year_divisions*(99.9-6.0)+1,2)];
    % apply prevalence in 5-year-olds to 0 to 6 year olds; apply prevalence in all ages to 6 to 99 year olds
    assert(isequal(size(StartPrev_byAgeGroups),[num_year_divisions*100 2]))

    NumSAg = StartPrev_byAgeGroups .* StartPop;
    NumNotSAg = (1 - StartPrev_byAgeGroups) .* StartPop;
elseif strcmp(seed_prev_string,'WHO')
    under_5_pos_vec_len = length(find(ages<=5.0));
    over_5_pos_vec_len = length(find(ages>5.0));
    assert(under_5_pos_vec_len+over_5_pos_vec_len==num_age_steps)
    StartPrev_byAgeGroups = [ ...
        repmat(demog.country_HBsAg_prevalences_by_ages_prevacc_young_old(1),under_5_pos_vec_len,2); ...
        repmat(demog.country_HBsAg_prevalences_by_ages_prevacc_young_old(2),over_5_pos_vec_len,2) ...
        ];
    assert(isequal(size(StartPrev_byAgeGroups),[(100*num_year_divisions) 2]))

    NumSAg = StartPrev_byAgeGroups .* StartPop;
    NumNotSAg = (1 - StartPrev_byAgeGroups) .* StartPop;
end    

X(1, :, :, 1) = NumNotSAg;
X(3, :, :, 1) = 0.5 * NumSAg;
X(4, :, :, 1) = 0.5 * NumSAg;

% Demography
% Prepare an index that will allow quick population of the mu vector from
% the demographic data input (uneven age-groupings)

MappingFromDataToParam = repmat(2:21,5*num_year_divisions,1);
MappingFromDataToParam = MappingFromDataToParam(:);
MappingFromDataToParam(1:num_year_divisions) = 1;
% MappingFromDataToParam gives the value in the mortality vectors (21 values
% corresponding to age groups 0--0, 1--4, 5--9, 10--14, ..., 80--84, 85--89, 90--94, 95--99) that should be
% used for the age groups 0, 0.1, 0.2, ..., 99.9

 
cov_InfantVacc_itt = demog.InfantVacc;
cov_BirthDose_itt = demog.BirthDose;
assert(all(cov_InfantVacc_itt >= 0))
assert(all(cov_InfantVacc_itt <= 1))
assert(all(cov_BirthDose_itt >= 0))
assert(all(cov_BirthDose_itt <= 1))
assert(isequal(size(cov_InfantVacc_itt),size(TimeSteps)))
assert(isequal(size(cov_BirthDose_itt),size(TimeSteps)))


assert(isequal(size(Prog),size(zeros(num_states, num_states)))); % Non-Age Specific Prog parameters stored as (from, to)



% Prepare for simulation

% prepare storage containers, for outputs once per year
% breakdowns by age/sex
[...
    Tot_Pop_1yr, Prev_Immune_Reactive_1yr, Prev_Chronic_Hep_B_1yr, Prev_Comp_Cirr_1yr, Prev_Decomp_Cirr_1yr, Prev_TDF_treat_1yr, NumSAg_1yr, NumSAg_chronic_1yr, yld_1yr, Prev_Deaths_1yr...
    ] = deal(-99 * ones(2, max(agegroups_1yr), T0 + 1));
[...
    Incid_chronic_all_1yr_approx,...
    Incid_Deaths_1yr_approx...
    ] = deal(-99 * ones(2, max(agegroups_1yr), T0 + 1));
							  
[...
    NewChronicCarriage, moving_btw_states, ...
    ] = deal(zeros(1, num_age_steps, 2, 2));

% single output per year
[Time, num_births_1yr] = deal(-99 * ones(1, T0 + 1));

% ----- Simulation -----

itt = 1; % itt increase every time i.e. every 0.1 years; goes from 1 to 2101 (length of TimeSteps)
OutputEventNum = 1; % OutputEventNum increase every year; goes from 1 to 212
moving_to_treatment = zeros(size(X));
initiated_treatment = false;
num_babies = 0;
babies_ChronicCarriage = 0;
female_multiplier = 0;
male_multiplier = 0;


for time = TimeSteps 
       
    
    % Update mortality and fertility rates
    mu = zeros(num_states, num_age_steps, 2, 2);
    mu(:, :, 1, :) = repmat(demog.MortalityRate_Women(OutputEventNum, MappingFromDataToParam), [num_states 1 2]);
    mu(:, :, 2, :) = repmat(demog.MortalityRate_Men(OutputEventNum, MappingFromDataToParam), [num_states 1 2]);
    % demog.MortalityRate_Women is a 212 x 21 matrix of mortality rates of 21 age groups (0--0, 1--4, 5--9, 10--14, ..., 90--94, 95--99) for every year from 1890 to 2101
    % OutputEventNum ranges from 1 to 212
    % agegroups_1yr = [1 1 1 ... 100 100 100], each number 10 times
    % selected vector copied across disease states and treatments
    % demog.fert is a 1000 x 212 matrix; ages in 0.1 year jumps versus 212 years
    fert = demog.fert(1:10:end, OutputEventNum);
    assert(isequal(size(fert),[100 1]))
    fert = repmat(fert',num_year_divisions,1);
    fert = fert(:);
    assert(isequal(size(fert),[100*num_year_divisions 1]))
    assert(length(demog.net_migration)==212)
    net_migration = demog.net_migration(OutputEventNum);
    sex_ratio = demog.sex_ratios(OutputEventNum);
    assert(length(net_migration)==1)
    net_migration = repmat(net_migration, [num_states num_age_steps 2 2]);


    
    % Compute Outputs once per year
    if rem(time, 1) == 0 % only saves variables in this "for" loop every 10 time steps (or once a year, since dt=0.1)
        Time(OutputEventNum) = time;
 
        
        % Rescale population sizes of each age group and gender
        if (time >= 1950) 
            ModelPop = squeeze(sum(sum(X([1:10 12:15],:,:,:), 1), 4));
            assert(isequal(size(ModelPop),[(100*num_year_divisions) 2]))
            % sum over disease state of alive people and treatment; ModelPop is 1000 x 2 i.e. age groups versus gender
            assert(isequal(size(demog.total_pop_female),[101 152]))
            % demog.total_pop_female is a 101 x 152 matrix of 152 years (1950 to 2101) for 101 age groups (0--0, 1--1, 2--2,..., 98--98, 99--99, 100--100)
            assert(isequal(size(demog.total_pop_male),[101 152]))
            col_index = time - 1949;
            MontaguPopFemale = demog.total_pop_female(1:100,col_index);
            % only want ages 0--99
            MontaguPopMale = demog.total_pop_male(1:100,col_index);
            MontaguPop = [MontaguPopFemale MontaguPopMale];
            assert(isequal(size(MontaguPop),[100 2]))
            MontaguPopExpand = MontaguPop(agegroups_1yr, :) * dt;
            % agegroups_1yr is a 1 x 1000 double; [1 1 ... 100 100], each number present 10 times
            % expanding MontaguPop from 1 year age steps to 0.1 year age steps
            % each age group repeated 10 times therefore divide each entry by 10
            assert(isequal(size(MontaguPopExpand),[(100*num_year_divisions) 2]))
            ScalerMat = MontaguPopExpand ./ ModelPop;
            ScalerMat(isnan(ScalerMat)) = 0;
            ScalerMat(isinf(ScalerMat)) = 0;
            pop_scaler = repmat(reshape(ScalerMat, [1 (100*num_year_divisions) 2]), [num_states 1 1 2]);
            % MontaguPopExpand is sizes of the current year's population over 0.1 year age steps; a 1000 x 2 matrix of ages versus gender
            % add an extra dimension and duplicate it for each disease state and treatment method
            X = X .* pop_scaler;
            % scale all parts of X, including dead people i.e. State 11
        end

    
        for k = 1:2 % genders

            
            for ag = 1:max(agegroups_1yr) % 1:100
 
                if OutputEventNum > 1
                
                    state_prev_vec = squeeze(sum(sum(X(:, agegroups_1yr == ag, k, :), 2), 4)); % k is gender
                    assert(isequal(size(state_prev_vec),[num_states 1]))

                    Tot_Pop_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec([1:10 12:15]));
                
                    Prev_Immune_Reactive_1yr(k, ag, OutputEventNum-1) = state_prev_vec(3);
                
                    Prev_Chronic_Hep_B_1yr(k, ag, OutputEventNum-1) = state_prev_vec(5);

                    Prev_Comp_Cirr_1yr(k, ag, OutputEventNum-1) = state_prev_vec(6);

                    Prev_Decomp_Cirr_1yr(k, ag, OutputEventNum-1) = state_prev_vec(7);
                
                    Prev_TDF_treat_1yr(k, ag, OutputEventNum-1) = state_prev_vec(10);

                    Prev_Deaths_1yr(k, ag, OutputEventNum-1) = state_prev_vec(11);

                    NumSAg_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec([2:8 10 12:15]));
                
                    NumSAg_chronic_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec([2:8 10 12:13]));

                    yld_1yr(k, ag, OutputEventNum-1) = sum( state_prev_vec .* demog.dwvec' );

                    Incid_chronic_all_1yr_approx(k,ag,OutputEventNum-1) = sum(sum(NewChronicCarriage(1, agegroups_1yr == ag, k, :), 2), 4);
                
                    Incid_Deaths_1yr_approx(k, ag, OutputEventNum-1) = sum(state_prev_vec .* Prog(:, 11));

                end
                
            end % end agegroups_1yr for loop
				
            if OutputEventNum > 1

                assert(Incid_chronic_all_1yr_approx(k, 1, OutputEventNum-1)==0) % 0-year-olds cannot get horizontal chronic infection (see FOI)							
                Incid_chronic_all_1yr_approx(1, 1, OutputEventNum-1) = female_multiplier * babies_ChronicCarriage;
                Incid_chronic_all_1yr_approx(2, 1, OutputEventNum-1) = male_multiplier * babies_ChronicCarriage;

            end
                
             
        end % end genders for loop
        
        if OutputEventNum > 1
            
            assert(length(num_babies)==1)
            num_births_1yr(OutputEventNum-1) = num_babies;
            
        end
        
        num_babies = 0;
        
        % Update counter
        OutputEventNum = OutputEventNum + 1;
        % increases every year
    end % end "rem(time, 1) == 0" if statement
    
    
    % Horizontal Transmission (infection and Chronic Carriage)
    
    FOI = zeros(1, num_age_steps, 2, 2);
    
    % i: Transmission Between 1y-5y olds
    FOI(1, i1y:(i5y - 1), :, :) = ...
        beta_U5_SAg_itt(itt) * sum(sum(sum(sum(X([4:8 13], i1y:(i5y - 1), :, :))))) / sum(sum(sum(sum(X([1:10 12:15], i1y:(i5y - 1), :, :))))) ...
        + beta_U5_EAg_itt(itt) * sum(sum(sum(sum(X([2:3 14:15], i1y:(i5y - 1), :, :))))) / sum(sum(sum(sum(X([1:10 12:15], i1y:(i5y - 1), :, :)))));
    
    % ii: Transmission between 1-15 year olds
    FOI(1, i1y:(i15y - 1), :, :) = FOI(1, i1y:(i15y - 1), :, :) + ...
        beta_1to15_SAg_itt(itt) * sum(sum(sum(sum(X([4:8 13], i1y:(i15y - 1), :, :))))) / sum(sum(sum(sum(X([1:10 12:15], i1y:(i15y - 1), :, :))))) ...
        + beta_1to15_EAg_itt(itt) * sum(sum(sum(sum(X([2:3 14:15], i1y:(i15y - 1), :, :))))) / sum(sum(sum(sum(X([1:10 12:15], i1y:(i15y - 1), :, :)))));
    
    % iii: Transmission Between 5+ and Adults (Assuming equal risks for all persons 5y-100y)
    FOI(1, i5y:end, :, :) = FOI(1, i5y:end, :, :) + ...
        beta_5plus_SAg_itt(itt) * sum(sum(sum(sum(X([4:8 13], i5y:end, :, :))))) / sum(sum(sum(sum(X([1:10 12:15], i5y:end, :, :))))) ...
        + beta_5plus_EAg_itt(itt) * sum(sum(sum(sum(X([2:3 14:15], i5y:end, :, :))))) / sum(sum(sum(sum(X([1:10 12:15], i5y:end, :, :)))));
    
    
    % Disease Progression
    next_X = X;
    for tr = 1:length(Transactions.From)
        transaction_vals = Transactions.Values{tr};
        transaction_vals = transaction_vals(:);
        moving_btw_states(1, :, :, :) = X(Transactions.From(tr), :, :, :) .* Transactions.Values{tr};
        next_X(Transactions.From(tr), :, :, :) = next_X(Transactions.From(tr), :, :, :) + dt * ( -moving_btw_states ); % move people out of "from" state
        next_X(Transactions.To(tr), :, :, :) = next_X(Transactions.To(tr), :, :, :)  + dt * ( +moving_btw_states ); % move people into "to" state
    end % end Disease Progression for loop
    % multiply by dt because quantity is calculated every 0.1 years therefore needs to be divided by 10



       
    % (Time-dependent) Baseline Transition to TDF-Treatment
    assert(squeeze(sum(sum(sum(sum(X([12 13], :, :, 1),1),2),3),4))==0)
    % (Time-dependent) Baseline Transition to TDF-Treatment
    if (time >= treat_start_year)
    % 2016 must be the first year with nonzero treatment 
    % therefore start treating from 2015.9 onwards since prevalence is recorded at the top of the loop

        if ~initiated_treatment
            num_in_treatment = sum(sum(sum(sum(X(10, :, :, :),1),2),3),4);
            assert(num_in_treatment==0) % no one is in treatment
            prev_pop = sum(sum(sum(sum(X([2:8 10 12:15], :, :, :),1),2),3),4); 

            total_num_to_move_to_treat = treat_coverage_in_2016 * prev_pop;
            eligible_pop = sum(sum(sum(sum(X([3 5 6 7], :, :, :),1),2),3),4); 
            assert(total_num_to_move_to_treat<eligible_pop)
            scaling_num = total_num_to_move_to_treat / eligible_pop;
            next_X([3 5 6 7],:,:,:)=next_X([3 5 6 7],:,:,:) - X([3 5 6 7],:,:,:) * scaling_num;
            % Every compartment in the eligible-for-treatment states in next_X must have a number subtracted from it 
            % such that the total number subtracted from the eligible-for-treatment states is in_treatment_2016
            % i.e. in_treatment_2016 = sum(sum(sum(sum(X([3 5 6 7], :, :, :),1),2),3),4) * scaling_num = sum(sum(sum(sum(X([3 5 6 7], :, :, :) * scaling_num,1),2),3),4)
            % Hence, scaling_num scales each compartment in X([3 5 6 7], :, :, :) such that X([3 5 6 7],:,:,:) * scaling_num subtracts the same proportion of people from each compartment in each of the eligible-for-treatment states in order to subtract a total of in_treatment_2016 from the eligible-for-treatment states.
            next_X(10,:,:,:) = next_X(10,:,:,:) + sum(X([3 5 6 7],:,:,:) * scaling_num,1);

            num_in_treatment = sum(sum(sum(sum(next_X(10, :, :, :),1),2),3),4);
            eligible_pop = sum(sum(sum(sum(X([3 5:7 10], :, :, :),1),2),3),4); 
            assert(num_in_treatment/eligible_pop >= treat_coverage_in_2016)
            % treatment coverage amongst treatment-eligible people will be greater than treatment coverage amongst HBsAg+ people, except if treatment coverage is 0
            treat_coverage_2016 = num_in_treatment / eligible_pop;

            initiated_treatment = true;
        else
            assert(initiated_treatment) % ensure that, each time this code is encountered, treatment has already been initiated

            assert(demog.PriorTDFTreatRate>=0)
            moving_to_treatment([3 5 6 7], :, :, :) = X([3 5 6 7], :, :, :) .* demog.PriorTDFTreatRate;
            next_X([3 5 6 7], :, :, :) = next_X([3 5 6 7], :, :, :) + dt * ( -moving_to_treatment([3 5 6 7], :, :, :) );
            next_X(10, :, :, :) = next_X(10, :, :, :) + dt * ( +sum(moving_to_treatment, 1) );
            assert(max(moving_to_treatment(:))>=0)
        end
    end % end treatment if statement
    

    
    
    % Infection process
    NewInfections = X(1, :, :, :) .* FOI;
    % number of susceptibles times FOI i.e. number of new infections within population, excluding babies (since FOI is 0 for babies)
    % 1 x num_age_steps x 2 x 2 double i.e. 1 x 1000 x 2 x 2 double
    NewChronicCarriage = NewInfections .* p_ChronicCarriage;
    SevereAcute = NewInfections * theta;
    NonsevereAcute = NewInfections - SevereAcute;
    
    % Transitions dependent on a state that does not have a number and is therefore not in Prog or Transactions
    % multiply by dt, since FOI is an annual rate
    next_X(1, :, :, :) = next_X(1, :, :, :) + dt * ( -NewInfections );
    next_X(14, :, :, :) = next_X(14, :, :, :) + dt * ( +NonsevereAcute );
    next_X(15, :, :, :) = next_X(15, :, :, :) + dt * ( +SevereAcute );
    
    
    % Infanct vaccination occurring at exactly six months
    % do not multiply by dt, since one is vaccinating cov_InfantVacc_itt(itt)% of people in next_X(1, i6mo, :, :), 
    % after which this cohort ages and moves to the next age bin
    % if divides all babies born in a year into 10 groups and vaccinatates cov_InfantVacc_itt(itt)% of each group, 
    % then one will have vaccinated cov_InfantVacc_itt(itt)% of all babies born in that year
    transfer_to_vacc = cov_InfantVacc_itt(itt) * next_X(1, i6mo, :, :) * Efficacy_InfantVacc; % the 0.95 represent a take-type vaccine efficacy of 95%.
    next_X(1, i6mo, :, :) = next_X(1, i6mo, :, :) - transfer_to_vacc;
    next_X(9, i6mo, :, :) = next_X(9, i6mo, :, :) + transfer_to_vacc;
    
        
    
    % Natural Mortality
    % Do not apply background mortality to the HBV deaths state, since we want people in all countries to be treated as if they would have lived until 84 if they had not died of HBV. This is done outside of the model in the main script. 
    mu(11,:,:,:,:)=0.0;
    next_X = next_X + dt * ( -next_X .* mu );
    % in the "for time=TimeSteps", which runs 10 times per year, therefore divide effect of mu by 10
    net_migration(11,:,:,:,:)=0.0;
    next_X = next_X + dt * ( +next_X .* net_migration );
    % "+" because net_migration = (number of immigrants - number of emigrants) / population size
    
    % Update Stocks
    X = next_X;
    
    % now age everyone
    X(:, 2:num_age_steps, :, :) = X(:, 1:(num_age_steps - 1), :, :);
    X(:, 1, :, :) = 0; % set number of babies to 0 (babies will be born next)
    
    
    % fill-out with new births in this time-step:
    births_toNonInfectiousWomen = sum( fert' .* sum(sum(X([1 9], :, 1, :), 1), 4) ); % Susecptible, Immune
    births_toHbEAgWomen = sum(fert' .* sum(sum(X([2:3 14:15], :, 1, :), 1), 4)); % Immune Tolerant, Immune Reactive
    births_toHbSAgWomen = sum(fert' .* sum(sum(X([4:8 13], :, 1, :), 1), 4)); % All other stages (other infected women)
    births_toTrWomen = sum(fert' .* sum(sum(X([10 12], :, 1, :), 1), 4)); % Women on Treatment
    
    
    births_Total = births_toNonInfectiousWomen + births_toHbEAgWomen + births_toHbSAgWomen + births_toTrWomen;
    assert(length(births_Total)==1)
    num_babies = num_babies + dt * births_Total;
    

    babies_ChronicCarriage = p_ChronicCarriage(1, 1, 1, 1) * ( ... % a 1 x 1 double
        ...
        births_toHbSAgWomen * (1 - cov_BirthDose_itt(itt)) * p_VerticalTransmission_HbSAg_NoIntv ...
        + births_toHbSAgWomen * cov_BirthDose_itt(itt) * p_VerticalTransmission_HbSAg_BirthDoseVacc ...
        ...
        + births_toHbEAgWomen * (1 - cov_BirthDose_itt(itt)) * p_VerticalTransmission_HbEAg_NoIntv ...
        + births_toHbEAgWomen * cov_BirthDose_itt(itt) * p_VerticalTransmission_HbEAg_BirthDoseVacc ...
        ...
        + births_toTrWomen * (1 - cov_BirthDose_itt(itt)) * p_VerticalTransmission_Tr_NoIntv ...
        + births_toTrWomen * cov_BirthDose_itt(itt) * p_VerticalTransmission_Tr_BirthDoseVacc ...
        );
    
    babies_NotChronicCarriage = births_Total - babies_ChronicCarriage;
    assert(length(babies_ChronicCarriage)==1)
    assert(length(babies_NotChronicCarriage)==1)

    
    female_multiplier = 1 / (1 + sex_ratio);
    male_multiplier = sex_ratio / (1 + sex_ratio);
    % sex_ratio is number of male births per one female birth
    % 1 / (1 + sex_ratio) + sex_ratio / (1 + sex_ratio) = 1, hence total number of babies not changed
    % male births -> 0 => sex_ratio -> 0 => male_multiplier -> 0/1 = 0
    % male births -> infinity => sex_ratio -> infinity => male_multiplier -> 1
    % female births -> 0 => sex_ratio -> infinity => female_multiplier -> 0
    % female births -> infinity => sex_ratio -> 0 => female_multiplier -> 1
    
    X(1, 1, 1, 1) = female_multiplier * dt * babies_NotChronicCarriage;  % Suscpetible babies
    X(2, 1, 1, 1) = female_multiplier * dt * babies_ChronicCarriage;     % Babies with chronic carriage
    
    X(1, 1, 2, 1) = male_multiplier * dt * babies_NotChronicCarriage;  % Suscpetible babies
    X(2, 1, 2, 1) = male_multiplier * dt * babies_ChronicCarriage;     % Babies with chronic carriage
      
    % increment the pointer
    itt = itt + 1;
    % increases every 0.1 years
    
end % end "time = TimeSteps" for loop

output.Time = Time; % 1 x (T0 + 1)
output.Tot_Pop_1yr = Tot_Pop_1yr; % 2 x 100 x (T0 + 1)
output.num_births_1yr = num_births_1yr; % 1 x (T0 + 1)
output.Incid_chronic_all_1yr_approx = Incid_chronic_all_1yr_approx; % 2 x 100 x (T0 + 1)
output.Prev_Immune_Reactive_1yr = Prev_Immune_Reactive_1yr; % 2 x 100 x (T0 + 1)
output.Prev_Chronic_Hep_B_1yr = Prev_Chronic_Hep_B_1yr; % 2 x 100 x (T0 + 1)
output.Prev_Comp_Cirr_1yr = Prev_Comp_Cirr_1yr; % 2 x 100 x (T0 + 1)
output.Prev_Decomp_Cirr_1yr = Prev_Decomp_Cirr_1yr; % 2 x 100 x (T0 + 1)
output.Prev_TDF_treat_1yr = Prev_TDF_treat_1yr; % 2 x 100 x (T0 + 1)
output.NumSAg_1yr = NumSAg_1yr; % 2 x 100 x (T0 + 1)
output.NumSAg_chronic_1yr = NumSAg_chronic_1yr; % 2 x 100 x (T0 + 1)
output.yld_1yr = yld_1yr; % 2 x 100 x (T0 + 1)
output.Incid_Deaths_1yr_approx = Incid_Deaths_1yr_approx; % 2 x 100 x (T0 + 1)
output.Prev_Deaths_1yr = Prev_Deaths_1yr; % 2 x 100 x (T0 + 1)

end % end function HBVmodel_PPT

