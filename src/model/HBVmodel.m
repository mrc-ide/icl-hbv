function output = HBVmodel(...
    demog,intvparams,effparams,seed_prev_string,num_states,num_year_divisions,dt,ages,num_age_steps,...
    start_year_simulations,end_year_simulations,T0,start_year_vacc_coverage,first_year_output,last_year_output,num_years_output,...
    theta,ECofactor,treat_start_year,treat_coverage_in_2016,p_ChronicCarriage,Prog,Transactions)


%%% demog = parametes to do with background demography and fitted parts of
%%% the model.
%%%

%% Establish basic simulation parameters
agegroups_5yr = 1 + floor(ages / 5); % categorises the ages into age-groups of 5 year width; 1 x 1000 double; [1 1 ... 20 20], each number present 50 times
agegroups_1yr = 1 + floor(ages); % categorises the ages into age-groups of 1 year width; 1 x 1000 double; [1 1 ... 100 100], each number present 10 times
					
i6mo = find(ages >= 0.5, 1); % markers for key age boundaries
i1y = find(ages >= 1, 1);
i5y = find(ages >= 5, 1);
i15y = find(ages >= 15, 1);

TimeSteps = start_year_simulations:dt:end_year_simulations; % 1 x 2101 double; [1890 1890.1 1890.2 ... 2099.8 2099.9 2100 2100.1 ... 2100.8 2100.9 2101]



% ---- Intervention Parameters ----

% * INFANT VACCINATION *
InfantVaccIntvCov = 0;                        % If this values is less than "current cov." no change is made. (0=> no change)
TScaleup_InfantVaccineIntv_start = 2022;          % The time at which scale-up starts for each respective intervention.
TScaleup_InfantVaccineIntv_finish = 2024;

% ** PMTCT **

    % Birth dose
    BirthDoseIntvCov = 0;                         % If this values is less than "current cov." no change is made. (0=> no change)
    TScaleup_BirthDoseVaccineIntv_start = 2022; 
    TScaleup_BirthDoseVaccineIntv_finish = 2024; 

    % PAP in addition to Birth dose
    % (Fraction of those who have BD that also get PAP)
    cov_BirthDoseAndTDF_EAgHighVL = 0;   
    cov_BirthDoseAndTDF_SAgHighVL = 0;
    cov_BirthDoseAndTDF_EAgLowVL = 0;
    cov_BirthDoseAndTDF_SAgLowVL = 0;

    % PAP instead of BD
    % (Fraction of those who do not get BD that do get PAP)
    cov_TDFOnly_EAgHighVL = 0;
    cov_TDFOnly_SAgHighVL = 0;
    cov_TDFOnly_EAgLowVL = 0;
    cov_TDFOnly_SAgLowVL = 0;

    TScaleup_PAP = 2025;    % This scale-up parameter pertains to both types of PAP usage (instantaneously)

% ** COUNTERFACTUALS **
TurnOffExistingInfantVacc = 0;        % Set infant vacc that has occured already to zero.
TurnOffExistingBirthDose = 0;         % Set birth dose vacc that has occured already to zero.
TurnOffExistingTreatment = 0;         % Set treatment that has occured already to zero.



%-------------
% PMTCT Asssumptions
% IMPORT THE EFFPARAMS (Now mandatory to give a full set of parameters)
FracEPosHighVL = effparams.FracEPosHighVL;
FracSPosHighVL = effparams.FracSPosHighVL;
p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL = effparams.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL;
p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL = effparams.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL;
pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc = effparams.pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc;
pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP = effparams.pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP;
pr_VerticalTransmission_HbSAgLowVL_PAP = effparams.pr_VerticalTransmission_HbSAgLowVL_PAP;
pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc = effparams.pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc;
pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP = effparams.pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP;
pr_VerticalTransmission_HbSAgHighVL_PAP = effparams.pr_VerticalTransmission_HbSAgHighVL_PAP;
pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc = effparams.pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc;
pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP = effparams.pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP;
pr_VerticalTransmission_HbEAgLowVL_PAP = effparams.pr_VerticalTransmission_HbEAgLowVL_PAP;
pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc = effparams.pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc;
pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP = effparams.pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP;
pr_VerticalTransmission_HbEAgHighVL_PAP = effparams.pr_VerticalTransmission_HbEAgHighVL_PAP;


% Confirm that the correct effparams have been entered:
effparams_required_params = {'FracEPosHighVL',...
                            'FracSPosHighVL',...
                            'p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL',...
                            'p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL',...
                            'pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc',...
                            'pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP',...
                            'pr_VerticalTransmission_HbSAgLowVL_PAP',...
                            'pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc',...
                            'pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP',...
                            'pr_VerticalTransmission_HbSAgHighVL_PAP',...
                            'pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc',...
                            'pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP',...
                            'pr_VerticalTransmission_HbEAgLowVL_PAP',...
                            'pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc',...
                            'pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP',...
                            'pr_VerticalTransmission_HbEAgHighVL_PAP'};
                    
for i=1:length(effparams_required_params)
    assert(exist(effparams_required_params{i},'var')>0)
end


% Compute p_VerticalTransmission_HbSAgHighVL_NoIntv and p_VerticalTransmission_HbSAgLowVL_NoIntv
p_HbSAg_av = demog.p_VerticalTransmission_HbSAg_NoIntv; % From fitting (this average rate to be preserved)
p_VerticalTransmission_HbSAgLowVL_NoIntv = p_HbSAg_av / ( (1-FracSPosHighVL) + p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL*FracSPosHighVL ); 
p_VerticalTransmission_HbSAgHighVL_NoIntv = min(1,p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL * p_VerticalTransmission_HbSAgLowVL_NoIntv);
assert(abs(p_HbSAg_av - (p_VerticalTransmission_HbSAgLowVL_NoIntv*(1-FracSPosHighVL) + p_VerticalTransmission_HbSAgHighVL_NoIntv*FracSPosHighVL))<0.001)
assert(p_VerticalTransmission_HbSAgLowVL_NoIntv>=0)
assert(p_VerticalTransmission_HbSAgLowVL_NoIntv<=1.0)
assert(p_VerticalTransmission_HbSAgHighVL_NoIntv>=0)
assert(p_VerticalTransmission_HbSAgHighVL_NoIntv<=1.0)

p_HbEAg_av = 0.9; % From assumption in the version of the model used in fitting  (this average rate to be preserved) [NB. if this varied in different region then specifiy that here!)]
p_VerticalTransmission_HbEAgLowVL_NoIntv = p_HbEAg_av / ( (1-effparams.FracEPosHighVL) + p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL*FracEPosHighVL ); 
p_VerticalTransmission_HbEAgHighVL_NoIntv = min(1,effparams.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL * p_VerticalTransmission_HbEAgLowVL_NoIntv);
assert(abs(p_HbEAg_av - (p_VerticalTransmission_HbEAgLowVL_NoIntv*(1-FracEPosHighVL) + p_VerticalTransmission_HbEAgHighVL_NoIntv*FracEPosHighVL))<0.001)
assert(p_VerticalTransmission_HbEAgHighVL_NoIntv>=0)
assert(p_VerticalTransmission_HbEAgHighVL_NoIntv<=1.0)

% Compute the transmission probabilities using the input probability ratios
p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc = pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc * p_VerticalTransmission_HbSAgLowVL_NoIntv;
p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP = pr_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP * p_VerticalTransmission_HbSAgLowVL_NoIntv;
p_VerticalTransmission_HbSAgLowVL_PAP = pr_VerticalTransmission_HbSAgLowVL_PAP * p_VerticalTransmission_HbSAgLowVL_NoIntv;

p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc = pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc * p_VerticalTransmission_HbSAgHighVL_NoIntv;
p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP = pr_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP * p_VerticalTransmission_HbSAgHighVL_NoIntv;
p_VerticalTransmission_HbSAgHighVL_PAP = pr_VerticalTransmission_HbSAgHighVL_PAP * p_VerticalTransmission_HbSAgHighVL_NoIntv;

p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc = pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc * p_VerticalTransmission_HbEAgLowVL_NoIntv;
p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP = pr_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP * p_VerticalTransmission_HbEAgLowVL_NoIntv;
p_VerticalTransmission_HbEAgLowVL_PAP = pr_VerticalTransmission_HbEAgLowVL_PAP * p_VerticalTransmission_HbEAgLowVL_NoIntv;

p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc = pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc * p_VerticalTransmission_HbEAgHighVL_NoIntv;
p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP = pr_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP * p_VerticalTransmission_HbEAgHighVL_NoIntv;
p_VerticalTransmission_HbEAgHighVL_PAP = pr_VerticalTransmission_HbEAgHighVL_PAP * p_VerticalTransmission_HbEAgHighVL_NoIntv;





% ------------------------------------------------------------------------
% Check that all the transmission probabilities make sense:

%   Boundary conditions
    assert(p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL>=1.0)
    assert(p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL>=1.0)

    for i = [1, 2, 5:length(effparams_required_params)]
        assert(eval(effparams_required_params{i})>=0.0)
        assert(eval(effparams_required_params{i})<=1.0)
    end
    
%   Within each e/vl cat, confirm appropriate ordering
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbSAgLowVL_NoIntv,p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbSAgLowVL_NoIntv,p_VerticalTransmission_HbSAgLowVL_PAP))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc,p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP))

    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbSAgHighVL_NoIntv,p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbSAgHighVL_NoIntv,p_VerticalTransmission_HbSAgHighVL_PAP))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc,p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP))
    
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgLowVL_NoIntv,p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgLowVL_NoIntv,p_VerticalTransmission_HbEAgLowVL_PAP))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc,p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP))

    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgHighVL_NoIntv,p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgHighVL_NoIntv,p_VerticalTransmission_HbEAgHighVL_PAP))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc,p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP))
    
%   That the high VL cat is always more transmissive than the low VL cat
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbSAgHighVL_NoIntv, p_VerticalTransmission_HbSAgLowVL_NoIntv))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgHighVL_NoIntv, p_VerticalTransmission_HbEAgLowVL_NoIntv))

    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc, p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc, p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc))

    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP, p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP, p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP))

    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbSAgHighVL_PAP, p_VerticalTransmission_HbSAgLowVL_PAP))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgHighVL_PAP, p_VerticalTransmission_HbEAgLowVL_PAP))

%   That the 'e cat' is always the same or more transmissive the the s cat
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgHighVL_NoIntv, p_VerticalTransmission_HbSAgHighVL_NoIntv))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgLowVL_NoIntv, p_VerticalTransmission_HbSAgLowVL_NoIntv))

    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc, p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc, p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc))

    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP, p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP, p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP)  )  

    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgHighVL_PAP, p_VerticalTransmission_HbSAgHighVL_PAP))
    assert(safe_greater_or_equal_to(p_VerticalTransmission_HbEAgLowVL_PAP, p_VerticalTransmission_HbSAgLowVL_PAP))
% ------------------------------------------------------------------------



    
%-------------


% ----- Allow overwriting of intv parameters  ---------
if nargin>1
    try InfantVaccIntvCov = intvparams.InfantVaccIntvCov; end
    try BirthDoseIntvCov = intvparams.BirthDoseIntvCov; end    
    try TScaleup_InfantVaccineIntv_start = intvparams.TScaleup_InfantVaccineIntv_start; end
    try TScaleup_InfantVaccineIntv_finish = intvparams.TScaleup_InfantVaccineIntv_finish; end
    try TScaleup_BirthDoseVaccineIntv_start = intvparams.TScaleup_BirthDoseVaccineIntv_start; end
    try TScaleup_BirthDoseVaccineIntv_finish = intvparams.TScaleup_BirthDoseVaccineIntv_finish; end

    try cov_BirthDoseAndTDF_EAgHighVL = intvparams.cov_BirthDoseAndTDF_EAgHighVL; end
    try cov_BirthDoseAndTDF_SAgHighVL = intvparams.cov_BirthDoseAndTDF_SAgHighVL; end
    try cov_BirthDoseAndTDF_EAgLowVL = intvparams.cov_BirthDoseAndTDF_EAgLowVL; end
    try cov_BirthDoseAndTDF_SAgLowVL = intvparams.cov_BirthDoseAndTDF_SAgLowVL; end
    try cov_TDFOnly_EAgHighVL = intvparams.cov_TDFOnly_EAgHighVL; end
    try cov_TDFOnly_SAgHighVL = intvparams.cov_TDFOnly_SAgHighVL; end
    try cov_TDFOnly_EAgLowVL = intvparams.cov_TDFOnly_EAgLowVL; end
    try cov_TDFOnly_SAgLowVL = intvparams.cov_TDFOnly_SAgLowVL; end
    try TScaleup_PAP = intvparams.TScaleup_PAP; end
    
    try TurnOffExistingInfantVacc = intvparams.TurnOffExistingInfantVacc; end
    try TurnOffExistingBirthDose = intvparams.TurnOffExistingBirthDose; end
end
% ------------------------------------------------





% ---- Load Epidemiological Parameters from fitting procedure ----
beta_U5 = demog.beta_U5;                                      % Rate of horizontal transmission between susceptible and infected persons - UNDER FIVE
beta_1to15 = demog.beta_1to15;                                % Rate of horizontal transmission between susceptible and infected persons - All Ages
beta_5plus = demog.beta_5plus;
ReducInTransmission = demog.ReducInTransmission;              % Fractional reduction in transmission
YearReducInTransmission = demog.YearReducInTransmission;      % Turning point year for reduction
DurReducInTransmission = 15;                                  % Time taken to complete change
PriorTDFTreatRate = demog.PriorTDFTreatRate;                  % Treatment rate since 2016



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
% demog.Pop_byAgeGroups_1950 is a 100 x 2 double of 100 age groups and 2 sexes
% agegroups_1yr is a 1 x 1000 double; [1 1 ... 100 100], each number present 10 times
% start population size like 1950 population
% expanding demog.Pop_byAgeGroups_1950 from 1 year age steps to 0.1 year age steps
% each age group repeated 10 times therefore divide each entry by 10

X = zeros(num_states, num_age_steps, 2, 2);
% dimensions: disease states, age, sex, accessible to treatment

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
    StartPrev_byAgeGroups = [...
        StartPrev_byAgeGroups(1:end-1,:); ...
        repmat(StartPrev_byAgeGroups(end,:),3,1)...
        ];
    % demog.HBsAg_prevs_middle_year_1 is a 18 x 2 double of age group by sex
    % demog.HBsAg_prevs_middle_year_1 age groups: 0--4 5--9 10--14 15--19 20--24 25--29 30--34 35--39 40--44 45--49 50--54 55--59 60--64 65--69 70--74 75--79 80--84 85+
    assert(isequal(size(StartPrev_byAgeGroups),[20 2]))

    NumSAg = StartPrev_byAgeGroups(agegroups_5yr, :) .* StartPop;
    % agegroups_5yr is a 1 x 1000 double; [1 1 ... 20 20], each number present 50 times
    % expanding StartPrev_byAgeGroups from 5 year age steps to 0.1 year age steps
    NumNotSAg = (1 - StartPrev_byAgeGroups(agegroups_5yr, :)) .* StartPop;
elseif strcmp(seed_prev_string,'CDA')
    StartPrev_byAgeGroups = [...
        repmat(demog.country_HBsAg_prevalences_by_ages_mid_1_young_old(1),num_year_divisions*(5.9-0.0)+1,2); ...
        repmat(demog.country_HBsAg_prevalences_by_ages_mid_1_young_old(2),num_year_divisions*(99.9-6.0)+1,2)...
        ];
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
end % end StartPop if statement

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



if end_year_simulations>2025

    last_available_year = 2019.0;
    index_last_available_year = last_available_year - 1979;
    assert(length(demog.InfantVacc)>=index_last_available_year)
    assert(length(demog.BirthDose)>=index_last_available_year)
    last_year_run = end_year_simulations;

    cov_InfantVacc = (TurnOffExistingInfantVacc==0) * demog.InfantVacc(1:index_last_available_year); % coverage of vaccination from 1980 to 2019
    cov_BirthDose = (TurnOffExistingBirthDose==0) * demog.BirthDose(1:index_last_available_year); 
    % 2019 is the last year of vaccination available from WUENIC, so scale-up to 90% is programmed to occur between 2019 and 2021

    xvals_vec = [[start_year_vacc_coverage:(last_available_year-1)] last_available_year TScaleup_InfantVaccineIntv_start TScaleup_InfantVaccineIntv_finish last_year_run];
    yvals_vec = [cov_InfantVacc(1:(end-1)) repmat(cov_InfantVacc(end),[1 2]) repmat(max(cov_InfantVacc(end),InfantVaccIntvCov),[1 2])];
    scaleup_years_vec = [TScaleup_InfantVaccineIntv_start TScaleup_InfantVaccineIntv_finish];
    cov_InfantVacc_itt_alt = make_coverage_vector(xvals_vec,yvals_vec,scaleup_years_vec,TimeSteps);
    cov_InfantVacc_itt = min(1,cov_InfantVacc_itt_alt);

    xvals_vec = [[start_year_vacc_coverage:(last_available_year-1)] last_available_year TScaleup_BirthDoseVaccineIntv_start TScaleup_BirthDoseVaccineIntv_finish last_year_run];
    yvals_vec = [cov_BirthDose(1:(end-1)) repmat(cov_BirthDose(end),[1 2]) repmat(max(cov_BirthDose(end),BirthDoseIntvCov),[1 2])];
    scaleup_years_vec = [TScaleup_BirthDoseVaccineIntv_start TScaleup_BirthDoseVaccineIntv_finish];
    cov_BirthDose_itt_alt = make_coverage_vector(xvals_vec,yvals_vec,scaleup_years_vec,TimeSteps);
    cov_BirthDose_itt = min(1,cov_BirthDose_itt_alt);


    % Coverage of PAP among those with BD
    xvals_vec = [start_year_vacc_coverage TScaleup_PAP-1 TScaleup_PAP last_year_run];
    yvals_vec = [0 0 cov_BirthDoseAndTDF_EAgHighVL cov_BirthDoseAndTDF_EAgHighVL];
    scaleup_years_vec = [TScaleup_PAP-1 TScaleup_PAP];
    cov_BirthDoseAndTDF_EAgHighVL_itt = interp1(xvals_vec,yvals_vec,TimeSteps,'linear','extrap');
    cov_BirthDoseAndTDF_EAgHighVL_itt = min(1,cov_BirthDoseAndTDF_EAgHighVL_itt);
    xvals_vec = [start_year_vacc_coverage TScaleup_PAP-1 TScaleup_PAP last_year_run];
    yvals_vec = [0 0 cov_BirthDoseAndTDF_EAgLowVL cov_BirthDoseAndTDF_EAgLowVL];
    scaleup_years_vec = [TScaleup_PAP-1 TScaleup_PAP];
    cov_BirthDoseAndTDF_EAgLowVL_itt = interp1(xvals_vec,yvals_vec,TimeSteps,'linear','extrap');
    cov_BirthDoseAndTDF_EAgLowVL_itt = min(1,cov_BirthDoseAndTDF_EAgLowVL_itt);
    xvals_vec = [start_year_vacc_coverage TScaleup_PAP-1 TScaleup_PAP last_year_run];
    yvals_vec = [0 0 cov_BirthDoseAndTDF_SAgHighVL cov_BirthDoseAndTDF_SAgHighVL];
    scaleup_years_vec = [TScaleup_PAP-1 TScaleup_PAP];
    cov_BirthDoseAndTDF_SAgHighVL_itt = interp1(xvals_vec,yvals_vec,TimeSteps,'linear','extrap');
    cov_BirthDoseAndTDF_SAgHighVL_itt = min(1,cov_BirthDoseAndTDF_SAgHighVL_itt);
    xvals_vec = [start_year_vacc_coverage TScaleup_PAP-1 TScaleup_PAP last_year_run];
    yvals_vec = [0 0 cov_BirthDoseAndTDF_SAgLowVL cov_BirthDoseAndTDF_SAgLowVL];
    scaleup_years_vec = [TScaleup_PAP-1 TScaleup_PAP];
    cov_BirthDoseAndTDF_SAgLowVL_itt = interp1(xvals_vec,yvals_vec,TimeSteps,'linear','extrap');
    cov_BirthDoseAndTDF_SAgLowVL_itt = min(1,cov_BirthDoseAndTDF_SAgLowVL_itt);

    % Coverage of PAP among those not with BD
    xvals_vec = [start_year_vacc_coverage TScaleup_PAP-1 TScaleup_PAP last_year_run];
    yvals_vec = [0 0 cov_TDFOnly_EAgHighVL cov_TDFOnly_EAgHighVL];
    scaleup_years_vec = [TScaleup_PAP-1 TScaleup_PAP];
    cov_TDFOnly_EAgHighVL_itt = interp1(xvals_vec,yvals_vec,TimeSteps,'linear','extrap');
    cov_TDFOnly_EAgHighVL_itt = min(1,cov_TDFOnly_EAgHighVL_itt);
    xvals_vec = [start_year_vacc_coverage TScaleup_PAP-1 TScaleup_PAP last_year_run];
    yvals_vec = [0 0 cov_TDFOnly_EAgLowVL cov_TDFOnly_EAgLowVL];
    scaleup_years_vec = [TScaleup_PAP-1 TScaleup_PAP];
    cov_TDFOnly_EAgLowVL_itt = interp1(xvals_vec,yvals_vec,TimeSteps,'linear','extrap');
    cov_TDFOnly_EAgLowVL_itt = min(1,cov_TDFOnly_EAgLowVL_itt);
    xvals_vec = [start_year_vacc_coverage TScaleup_PAP-1 TScaleup_PAP last_year_run];
    yvals_vec = [0 0 cov_TDFOnly_SAgHighVL cov_TDFOnly_SAgHighVL];
    scaleup_years_vec = [TScaleup_PAP-1 TScaleup_PAP];
    cov_TDFOnly_SAgHighVL_itt = interp1(xvals_vec,yvals_vec,TimeSteps,'linear','extrap');
    cov_TDFOnly_SAgHighVL_itt = min(1,cov_TDFOnly_SAgHighVL_itt);
    xvals_vec = [start_year_vacc_coverage TScaleup_PAP-1 TScaleup_PAP last_year_run];
    yvals_vec = [0 0 cov_TDFOnly_SAgLowVL cov_TDFOnly_SAgLowVL];
    scaleup_years_vec = [TScaleup_PAP-1 TScaleup_PAP];
    cov_TDFOnly_SAgLowVL_itt = interp1(xvals_vec,yvals_vec,TimeSteps,'linear','extrap');
    cov_TDFOnly_SAgLowVL_itt = min(1,cov_TDFOnly_SAgLowVL_itt);

else

    assert(all(TimeSteps)<start_year_vacc_coverage) % vaccination does not start before 1980


    [cov_InfantVacc_itt, ...
        cov_BirthDose_itt] = ...
        deal(zeros(size(TimeSteps)));


    % Coverage of PAP among those with BD
    [cov_BirthDoseAndTDF_EAgHighVL_itt, ...
        cov_BirthDoseAndTDF_EAgLowVL_itt, ...
        cov_BirthDoseAndTDF_SAgHighVL_itt, ...
        cov_BirthDoseAndTDF_SAgLowVL_itt] = ...
        deal(zeros(size(TimeSteps)));

    % Coverage of PAP among those not with BD
    [cov_TDFOnly_EAgHighVL_itt, ...
        cov_TDFOnly_EAgLowVL_itt, ...
        cov_TDFOnly_SAgHighVL_itt, ...
        cov_TDFOnly_SAgLowVL_itt] = ...
        deal(zeros(size(TimeSteps)));

end % end_year_simulations>2025 if statement

assert(all(cov_InfantVacc_itt >= 0) && all(cov_InfantVacc_itt <= 1))
assert(all(cov_BirthDose_itt >= 0) && all(cov_BirthDose_itt <= 1))
assert(isequal(size(cov_InfantVacc_itt),size(TimeSteps)) && isequal(size(cov_BirthDose_itt),size(TimeSteps)))


assert(isequal(size(Prog),[num_states, num_states])); % Non-Age Specific Prog parameters stored as (from, to)



% Prepare for simulation

% prepare storage containers, for outputs once per year
% breakdowns by age/sex

[NumSAg_5yr, PrevEAg_of_SAg_5yr] = deal(-99 * ones(2, max(agegroups_5yr), T0 + 1)); % 1890 + T0 gives T0 + 1 years in total

Incid_chronic_all_5yr_approx = zeros(2, max(agegroups_5yr), T0 + 1);

[...
    Tot_Pop_1yr, Prev_Decomp_Cirr_1yr, Prev_Liver_Cancer_1yr, NumSAg_1yr, NumEAg_chronic_1yr, NumEAg_chronic_acute_1yr, yld_1yr, Prev_Deaths_1yr...
    ] = deal(-99 * ones(2, max(agegroups_1yr), T0 + 1));

Incid_Deaths_1yr_approx = deal(-99 * ones(2, max(agegroups_1yr), T0 + 1));
							  
% single output per year

[...
    Time, RateInfantVacc, RateBirthDoseVacc, RatePeripartumTreatment, NumDecompCirr, NumLiverCancer, PregnantWomenNeedToScreen, HBVPregnantWomenNeedToEvaluate...
    ] = deal(-99 * ones(1, T0 + 1));

[num_births_toHbEAgWomenHVL_1yr_approx, num_births_toHbEAgWomenLVL_1yr_approx, num_births_toHbSAgWomenHVL_1yr_approx, num_births_toHbSAgWomenLVL_1yr_approx, ... 
    num_births_1yr_approx, ...
    num_births_chronic_HbEAgWomenHVL_1yr_approx, num_births_chronic_HbEAgWomenLVL_1yr_approx, num_births_chronic_HbSAgWomenHVL_1yr_approx, num_births_chronic_HbSAgWomenLVL_1yr_approx, ... 
    Incid_babies_chronic_1yr_approx, ...
    PeripartumTreatment_HbEAg_HighVL_approx, PeripartumTreatment_HbEAg_LowVL_approx, PeripartumTreatment_HbSAg_HighVL_approx, PeripartumTreatment_HbSAg_LowVL_approx...
    ] = deal(-99 * ones(1, T0 + 1));

% ----- Simulation -----

itt = 1; % itt increase every time i.e. every 0.1 years; goes from 1 to 2101 (length of TimeSteps)
OutputEventNum = 1; % OutputEventNum increase every year; goes from 1 to 212
moving_to_treatment = zeros(size(X));
initiated_treatment = false;
moving_btw_states = zeros(1, num_age_steps, 2, 2);
into_babies_chronic = zeros(2,1);
transfer_to_vacc = zeros(1, 1, 2, 2);
[FOI, NewChronicCarriage] = deal(zeros(1, num_age_steps, 2, 2));

[births_toHbEAgWomenHighVL, births_toHbEAgWomenLowVL, births_toHbSAgWomenHighVL, births_toHbSAgWomenLowVL,...
    births_Total,...
    babiesChronic_from_HbEAgWomenHighVL, babiesChronic_from_HbEAgWomenLowVL, babiesChronic_from_HbSAgWomenHighVL, babiesChronic_from_HbSAgWomenLowVL,...
    ratebirthdoses,...
    babies_ChronicCarriage,...
    pregnantWomenNeedToScreen,...
    num_mothers_PAP_HbEAg_HighVL, num_mothers_PAP_HbEAg_LowVL, num_mothers_PAP_HbSAg_HighVL, num_mothers_PAP_HbSAg_LowVL,...
    RateOfPAPInitiation,...
    HBVPositivePregnantWomenAtANC,...
    female_multiplier, male_multiplier] = deal(0);




for time = TimeSteps 
       
    
    % Update mortality and fertility rates
    mu = zeros(num_states, num_age_steps, 2, 2);
    assert(isequal(size(demog.MortalityRate_Women),size(demog.MortalityRate_Men)))
    mu(:, :, 1, :) = repmat(demog.MortalityRate_Women(OutputEventNum, MappingFromDataToParam), [num_states 1 2]);
    mu(:, :, 2, :) = repmat(demog.MortalityRate_Men(OutputEventNum, MappingFromDataToParam), [num_states 1 2]);
    % demog.MortalityRate_Women is a 212 x 21 matrix of mortality rates of 21 age groups (0--0, 1--4, 5--9, 10--14, ..., 90--94, 95--99) for every year from 1890 to 2101
    % OutputEventNum ranges from 1 to 212
    % agegroups_1yr = [1 1 1 ... 100 100 100], each number 10 times
    % selected vector copied across disease states and treatments
    assert(isequal(size(demog.fert),[(100*num_year_divisions) 212]))
    % demog.fert is a 1000 x 212 matrix; ages in 0.1 year jumps versus 212 years
    fert = demog.fert(1:10:end, OutputEventNum);
    assert(isequal(size(fert),[100 1]))
    fert = repmat(fert',num_year_divisions,1);
    fert = fert(:);
    assert(isequal(size(fert),[100*num_year_divisions 1]))
    assert(length(demog.net_migration)==212)
    net_migration = demog.net_migration(OutputEventNum);
    assert(length(net_migration)==1)
    net_migration = repmat(net_migration, [num_states num_age_steps 2 2]);
    sex_ratio = demog.sex_ratios(OutputEventNum);

    
    % Compute Outputs once per year
    if rem(time, 1) == 0 % only saves variables in this "for" loop every 10 time steps (or once a year, since dt=0.1)


        Time(OutputEventNum) = time;
 
        
        % Rescale population sizes of each age group and sex
        if (time >= 1950) 
            ModelPop = squeeze(sum(sum(X([1:10 12:15],:,:,:), 1), 4));
            assert(isequal(size(ModelPop),[(100*num_year_divisions) 2]))
            % sum over disease state of alive people and treatment; ModelPop is 1000 x 2 i.e. age groups versus sex
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
            % MontaguPopExpand is sizes of the current year's population over 0.1 year age steps; a 1000 x 2 matrix of ages versus sex
            % add an extra dimension and duplicate it for each disease state and treatment method
            X = X .* pop_scaler;
            % scale all parts of X, including dead people i.e. State 11
        end
        
        
        for k = 1:2 % sexes
                
            for ag = 1:max(agegroups_5yr) % 1:20
														
                if OutputEventNum > 1
                    % model results assigned to a particular year at the beginning of that year, after which they are zeroed

                    NumSAg_5yr(k, ag, OutputEventNum-1) = sum(sum(sum(X([2:8 10 12:15], agegroups_5yr == ag, k, :))));
                
                    if NumSAg_5yr(k, ag, OutputEventNum-1)>0
                        PrevEAg_of_SAg_5yr(k,ag,OutputEventNum-1) = sum(sum(sum(X(2:3, agegroups_5yr == ag, k, :)))) / NumSAg_5yr(k, ag, OutputEventNum-1);
                        % Note that this is prevalence of e+ among s+
                    else
                        assert(NumSAg_5yr(k, ag, OutputEventNum-1)==0)
                        assert(sum(sum(sum(X(2:3, agegroups_5yr == ag, k, :))))==0)
                        PrevEAg_of_SAg_5yr(k,ag,OutputEventNum-1) = 0;
                    end

                    Incid_chronic_all_5yr_approx(k, ag, OutputEventNum-1) = sum(sum(NewChronicCarriage(1, agegroups_5yr == ag, k, :), 2), 4);

                end
                
            end

            
            for ag = 1:max(agegroups_1yr) % 1:100
 
                if OutputEventNum > 1
                
                    state_prev_vec = squeeze(sum(sum(X(:, agegroups_1yr == ag, k, :), 2), 4)); % k is sex
                    assert(isequal(size(state_prev_vec),[num_states 1]))

                    Tot_Pop_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec([1:10 12:15]));

                    Prev_Decomp_Cirr_1yr(k, ag, OutputEventNum-1) = state_prev_vec(7);

                    Prev_Liver_Cancer_1yr(k, ag, OutputEventNum-1) = state_prev_vec(8);

                    Prev_Deaths_1yr(k, ag, OutputEventNum-1) = state_prev_vec(11);

                    NumSAg_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec([2:8 10 12:15]));

                    NumEAg_chronic_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec([2:3]));

                    NumEAg_chronic_acute_1yr(k, ag, OutputEventNum-1) = sum(state_prev_vec([2:3 14:15]));

                    yld_1yr(k, ag, OutputEventNum-1) = sum( state_prev_vec .* demog.dwvec' );
                
                    Incid_Deaths_1yr_approx(k, ag, OutputEventNum-1) = sum(state_prev_vec .* Prog(:, 11));

                end
                
            end % end agegroups_1yr for loop
				
            if OutputEventNum > 1

                Incid_chronic_all_5yr_approx(1, 1, OutputEventNum-1) = female_multiplier * babies_ChronicCarriage;
                Incid_chronic_all_5yr_approx(2, 1, OutputEventNum-1) = male_multiplier * babies_ChronicCarriage;

            end
                
             
        end % end k for loop
        
        if OutputEventNum > 1
                
            state_prev_vec = squeeze(sum(sum(sum(X,2),3),4)); % 15 x 1
        
            NumDecompCirr(OutputEventNum-1) = state_prev_vec(7);
        
            NumLiverCancer(OutputEventNum-1) = state_prev_vec(8);

            num_births_toHbEAgWomenHVL_1yr_approx(OutputEventNum-1) = births_toHbEAgWomenHighVL;

            num_births_toHbEAgWomenLVL_1yr_approx(OutputEventNum-1) = births_toHbEAgWomenLowVL;

            num_births_toHbSAgWomenHVL_1yr_approx(OutputEventNum-1) = births_toHbSAgWomenHighVL;

            num_births_toHbSAgWomenLVL_1yr_approx(OutputEventNum-1) = births_toHbSAgWomenLowVL;

            num_births_1yr_approx(OutputEventNum-1) = births_Total;

            num_births_chronic_HbEAgWomenHVL_1yr_approx(OutputEventNum-1) = babiesChronic_from_HbEAgWomenHighVL;

            num_births_chronic_HbEAgWomenLVL_1yr_approx(OutputEventNum-1) = babiesChronic_from_HbEAgWomenLowVL;

            num_births_chronic_HbSAgWomenHVL_1yr_approx(OutputEventNum-1) = babiesChronic_from_HbSAgWomenHighVL;

            num_births_chronic_HbSAgWomenLVL_1yr_approx(OutputEventNum-1) = babiesChronic_from_HbSAgWomenLowVL;

            Incid_babies_chronic_1yr_approx(OutputEventNum-1) = babies_ChronicCarriage;
            
            RateBirthDoseVacc(OutputEventNum-1) = ratebirthdoses;
        
            RateInfantVacc(OutputEventNum-1) = squeeze(sum(sum(transfer_to_vacc,3),4)) * num_year_divisions;
             
            PregnantWomenNeedToScreen(OutputEventNum-1) = pregnantWomenNeedToScreen;
        
            PeripartumTreatment_HbEAg_HighVL_approx(OutputEventNum-1) = num_mothers_PAP_HbEAg_HighVL;

            PeripartumTreatment_HbEAg_LowVL_approx(OutputEventNum-1) = num_mothers_PAP_HbEAg_LowVL;

            PeripartumTreatment_HbSAg_HighVL_approx(OutputEventNum-1) = num_mothers_PAP_HbSAg_HighVL;

            PeripartumTreatment_HbSAg_LowVL_approx(OutputEventNum-1) = num_mothers_PAP_HbSAg_LowVL;

            RatePeripartumTreatment(OutputEventNum-1) = RateOfPAPInitiation;
        
            HBVPregnantWomenNeedToEvaluate(OutputEventNum-1) = HBVPositivePregnantWomenAtANC;
      
        end
        
        into_babies_chronic = zeros(2,1);    
        
        
        % Update counter
        OutputEventNum = OutputEventNum + 1; % increases every year


    end % end "rem(time, 1) == 0" if statement
    
    
    
    % Horizontal Transmission (infection and Chronic Carriage)
    
    FOI = zeros(1, num_age_steps, 2, 2);
    
    % i: Transmission Between Children (1y-5y olds)
    FOI(1, i1y:(i5y - 1), :, :) = ...
        beta_U5_SAg_itt(itt) * sum(sum(sum(sum(X([4:8 13], i1y:(i5y - 1), :, :))))) / sum(sum(sum(sum(X([1:10 12:15], i1y:(i5y - 1), :, :))))) ...
        + beta_U5_EAg_itt(itt) * sum(sum(sum(sum(X([2:3 14:15], i1y:(i5y - 1), :, :))))) / sum(sum(sum(sum(X([1:10 12:15], i1y:(i5y - 1), :, :)))));
    
    % ii. Transmission between 1-15 year olds
    FOI(1, i1y:(i15y - 1), :, :) = FOI(1, i1y:(i15y - 1), :, :) + ...
        beta_1to15_SAg_itt(itt) * sum(sum(sum(sum(X([4:8 13], i1y:(i15y - 1), :, :))))) / sum(sum(sum(sum(X([1:10 12:15], i1y:(i15y - 1), :, :))))) ...
        + beta_1to15_EAg_itt(itt) * sum(sum(sum(sum(X([2:3 14:15], i1y:(i15y - 1), :, :))))) / sum(sum(sum(sum(X([1:10 12:15], i1y:(i15y - 1), :, :)))));
    
    % iii.: Transmission Between 5+ and Adults (Assuming equal risks for all persons 5y-100y)
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
       


    assert(squeeze(sum(sum(sum(sum(X([12 13], :, :, 1),1),2),3),4))==0)
    % (Time-dependent) Baseline Transition to TDF-Treatment
    if (time >= treat_start_year && treat_coverage_in_2016 > 0)
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
    assert(isequal(size(X(1, :, :, :)),size(FOI)))
    NewInfections = X(1, :, :, :) .* FOI;
    % number of susceptibles times FOI i.e. number of new infections within population, excluding babies (since FOI is 0 for babies)
    % 1 x num_age_steps x 2 x 2 double i.e. 1 x 1000 x 2 x 2 double
    assert(isequal(size(NewInfections),size(p_ChronicCarriage)))
    NewChronicCarriage = NewInfections .* p_ChronicCarriage;
    SevereAcute = NewInfections * theta;
    NonsevereAcute = NewInfections - SevereAcute;
    
    % Transitions dependent on a state that does not have a number and is therefore not in Prog or Transactions
    % multiply by dt, since FOI is an annual rate
    next_X(1, :, :, :) = next_X(1, :, :, :) + dt * ( -NewInfections );
    next_X(14, :, :, :) = next_X(14, :, :, :) + dt * ( +NonsevereAcute );
    next_X(15, :, :, :) = next_X(15, :, :, :) + dt * ( +SevereAcute );
    
    
    % Infant vaccination occurring at exactly six months
    % do not multiply by dt, since one is vaccinating cov_InfantVacc_itt(itt)% of people in next_X(1, i6mo, :, :), 
    % after which this cohort ages and moves to the next age bin
    % if divides all babies born in a year into 10 groups and vaccinatates cov_InfantVacc_itt(itt)% of each group, 
    % then one will have vaccinated cov_InfantVacc_itt(itt)% of all babies born in that year
    transfer_to_vacc = cov_InfantVacc_itt(itt) * next_X(1, i6mo, :, :) * demog.Efficacy_InfantVacc;
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
    
    
    
    % Births and MTCT:
    
    % fill-out with new births at the rate they occur
    births_toNonInfectiousWomen = sum( fert' .* sum(sum(X([1 9 10 12],:,1,:),1),4) ); %Susecptible, Immune, or on Treatment
    births_toHbEAgWomenHighVL = FracEPosHighVL * sum(fert' .* sum(sum(X([2 3 14 15],:,1,:),1),4)); %Immune Tolerant, Immune Reactive
    births_toHbEAgWomenLowVL = (1-FracEPosHighVL) * sum(fert' .* sum(sum(X([2 3 14 15],:,1,:),1),4)); %Immune Tolerant, Immune Reactive
    births_toHbSAgWomenHighVL = FracSPosHighVL * sum(fert' .* sum(sum(X([4 5 6 7 8 13],:,1,:),1),4)); %All other stages (other infected women)
    births_toHbSAgWomenLowVL = (1-FracSPosHighVL) * sum(fert' .* sum(sum(X([4 5 6 7 8 13],:,1,:),1),4)); %All other stages (other infected women)
    births_Total = births_toNonInfectiousWomen + births_toHbEAgWomenHighVL + births_toHbEAgWomenLowVL + births_toHbSAgWomenHighVL + births_toHbSAgWomenLowVL;

    babiesChronic_from_HbEAgWomenHighVL = p_ChronicCarriage(1, 1, 1, 1) * births_toHbEAgWomenHighVL * (...
                (1-cov_BirthDose_itt(itt)) * (1-cov_TDFOnly_EAgHighVL_itt(itt)) * p_VerticalTransmission_HbEAgHighVL_NoIntv  ...
            +   (1-cov_BirthDose_itt(itt)) * cov_TDFOnly_EAgHighVL_itt(itt) * p_VerticalTransmission_HbEAgHighVL_PAP  ...
            +   cov_BirthDose_itt(itt) * (1-cov_BirthDoseAndTDF_EAgHighVL_itt(itt)) * p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc  ...
            +   cov_BirthDose_itt(itt) * cov_BirthDoseAndTDF_EAgHighVL_itt(itt) * p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP  ...
            );
    % number of chronic babies born to HVL HBeAg+ pregnant women =
    %     (probability of infection becoming chronic in babies) * (number of births to HVL HBeAg+ pregnant women) *
    %     (
    %       probability of a baby that does not receive BD born to a HVL HBeAg+ pregnant woman who does not receive PAP being infected
    %     + probability of a baby that does not receive BD born to a HVL HBeAg+ pregnant woman who receives PAP being infected
    %     + probability of a baby that receives BD born to a HVL HBeAg+ pregnant woman who does not receive PAP being infected
    %     + probability of a baby that receives BD born to a HVL HBeAg+ pregnant woman who receives PAP being infected
    %     )
                                                                       
    babiesChronic_from_HbEAgWomenLowVL = p_ChronicCarriage(1, 1, 1, 1) * births_toHbEAgWomenLowVL * (...
                (1-cov_BirthDose_itt(itt)) * (1-cov_TDFOnly_EAgLowVL_itt(itt)) * p_VerticalTransmission_HbEAgLowVL_NoIntv  ...
            +   (1-cov_BirthDose_itt(itt)) * cov_TDFOnly_EAgLowVL_itt(itt) * p_VerticalTransmission_HbEAgLowVL_PAP  ...
            +   cov_BirthDose_itt(itt) * (1-cov_BirthDoseAndTDF_EAgLowVL_itt(itt)) * p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc  ...
            +   cov_BirthDose_itt(itt) * cov_BirthDoseAndTDF_EAgLowVL_itt(itt) * p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP  ...
            );

    babiesChronic_from_HbSAgWomenHighVL = p_ChronicCarriage(1, 1, 1, 1) * births_toHbSAgWomenHighVL * (...
                (1-cov_BirthDose_itt(itt)) * (1-cov_TDFOnly_SAgHighVL_itt(itt)) * p_VerticalTransmission_HbSAgHighVL_NoIntv  ...
            +   (1-cov_BirthDose_itt(itt)) * cov_TDFOnly_SAgHighVL_itt(itt) * p_VerticalTransmission_HbSAgHighVL_PAP  ...
            +   cov_BirthDose_itt(itt) * (1-cov_BirthDoseAndTDF_SAgHighVL_itt(itt)) * p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc  ...
            +   cov_BirthDose_itt(itt) * cov_BirthDoseAndTDF_SAgHighVL_itt(itt) * p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP  ...
            );
                                                                       
    babiesChronic_from_HbSAgWomenLowVL = p_ChronicCarriage(1, 1, 1, 1) * births_toHbSAgWomenLowVL * (...
                (1-cov_BirthDose_itt(itt)) * (1-cov_TDFOnly_SAgLowVL_itt(itt)) * p_VerticalTransmission_HbSAgLowVL_NoIntv  ...
            +   (1-cov_BirthDose_itt(itt)) * cov_TDFOnly_SAgLowVL_itt(itt) * p_VerticalTransmission_HbSAgLowVL_PAP  ...
            +   cov_BirthDose_itt(itt) * (1-cov_BirthDoseAndTDF_SAgLowVL_itt(itt)) * p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc  ...
            +   cov_BirthDose_itt(itt) * cov_BirthDoseAndTDF_SAgLowVL_itt(itt) * p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP  ...
            );
                                                                       
                                                                     
    num_mothers_PAP_HbEAg_HighVL = births_toHbEAgWomenHighVL * ( (1-cov_BirthDose_itt(itt))*cov_TDFOnly_EAgHighVL_itt(itt) + cov_BirthDose_itt(itt)*cov_BirthDoseAndTDF_EAgHighVL_itt(itt) );
    % HVL HbEAg+ pregnant women who get PAP = (HVL HbEAg+ pregnant women who get PAP and whose babies do not receive BD) + (HVL HbEAg+ pregnant women who get PAP and whose babies receive BD),
    % where number of births is used to approximate number of mothers (a woman can have twins, which makes number of births not equal to number of mothers)
    num_mothers_PAP_HbEAg_LowVL = births_toHbEAgWomenLowVL  * ( (1-cov_BirthDose_itt(itt))*cov_TDFOnly_EAgLowVL_itt(itt) + cov_BirthDose_itt(itt)*cov_BirthDoseAndTDF_EAgLowVL_itt(itt)  );
    % LVL HbEAg+ pregnant women who get PAP, where number of births is used to approximate number of mothers
    num_mothers_PAP_HbSAg_HighVL = births_toHbSAgWomenHighVL * ( (1-cov_BirthDose_itt(itt))*cov_TDFOnly_SAgHighVL_itt(itt) + cov_BirthDose_itt(itt)*cov_BirthDoseAndTDF_SAgHighVL_itt(itt)  );
    % HVL HbEAg- pregnant women who get PAP, where number of births is used to approximate number of mothers
    num_mothers_PAP_HbSAg_LowVL = births_toHbSAgWomenLowVL  * ( (1-cov_BirthDose_itt(itt))*cov_TDFOnly_SAgLowVL_itt(itt) + cov_BirthDose_itt(itt)*cov_BirthDoseAndTDF_SAgLowVL_itt(itt)  );
    % LVL HbEAg- pregnant women who get PAP, where number of births is used to approximate number of mothers
    RateOfPAPInitiation = num_mothers_PAP_HbEAg_HighVL + num_mothers_PAP_HbEAg_LowVL + num_mothers_PAP_HbSAg_HighVL + num_mothers_PAP_HbSAg_LowVL;
    % pregnant women who get PAP

    ratebirthdoses = births_Total * cov_BirthDose_itt(itt);

    babies_ChronicCarriage = babiesChronic_from_HbEAgWomenHighVL + babiesChronic_from_HbEAgWomenLowVL + babiesChronic_from_HbSAgWomenHighVL + babiesChronic_from_HbSAgWomenLowVL;

    babies_NotChronicCarriage = births_Total - babies_ChronicCarriage;

    tmp_MTCTRate_SPosPregWomen = (babiesChronic_from_HbSAgWomenHighVL + babiesChronic_from_HbSAgWomenLowVL) / (births_toHbSAgWomenLowVL + births_toHbSAgWomenHighVL);

    tmp_MTCTRate_EPosPregWomen = (babiesChronic_from_HbEAgWomenHighVL + babiesChronic_from_HbEAgWomenLowVL) / (births_toHbEAgWomenLowVL + births_toHbEAgWomenHighVL);

    tmp_MTCTRate_AllPosPregWomen = ...
                    (babiesChronic_from_HbEAgWomenHighVL + babiesChronic_from_HbEAgWomenLowVL + babiesChronic_from_HbSAgWomenHighVL + babiesChronic_from_HbSAgWomenLowVL) / ...
                    (births_toHbSAgWomenLowVL + births_toHbSAgWomenHighVL + births_toHbEAgWomenLowVL + births_toHbEAgWomenHighVL);

    pregnantWomenNeedToScreen = births_Total * max([cov_BirthDoseAndTDF_EAgHighVL_itt(itt),...
                                                    cov_BirthDoseAndTDF_EAgLowVL_itt(itt),...
                                                    cov_BirthDoseAndTDF_SAgHighVL_itt(itt),...
                                                    cov_BirthDoseAndTDF_SAgLowVL_itt(itt), ...
                                                    cov_TDFOnly_EAgHighVL_itt(itt),...
                                                    cov_TDFOnly_EAgLowVL_itt(itt),...
                                                    cov_TDFOnly_SAgHighVL_itt(itt),...
                                                    cov_TDFOnly_SAgLowVL_itt(itt)]);


    HBVPositivePregnantWomenAtANC = (births_toHbSAgWomenHighVL + births_toHbSAgWomenLowVL + births_toHbEAgWomenHighVL + births_toHbEAgWomenLowVL) ...        %added 13/8/19
                                    * max([cov_BirthDoseAndTDF_EAgHighVL_itt(itt),...
                                           cov_BirthDoseAndTDF_EAgLowVL_itt(itt),...
                                           cov_BirthDoseAndTDF_SAgHighVL_itt(itt),...
                                           cov_BirthDoseAndTDF_SAgLowVL_itt(itt),...
                                           cov_TDFOnly_EAgHighVL_itt(itt),...
                                           cov_TDFOnly_EAgLowVL_itt(itt),...
                                           cov_TDFOnly_SAgHighVL_itt(itt),...
                                           cov_TDFOnly_SAgLowVL_itt(itt)]);

    
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
  
    assert(isequal(size(into_babies_chronic),[2 1]))
    assert(isequal(size(squeeze(X(2, 1, :, 1))),[2 1]))
    into_babies_chronic = into_babies_chronic + squeeze(X(2, 1, :, 1));

  
    % Prevalence HBsAg among pregnant women
    tmp_PrevalenceAmongPregnantWomen = ...
        (births_toHbSAgWomenLowVL + births_toHbSAgWomenHighVL + births_toHbEAgWomenLowVL + births_toHbEAgWomenHighVL) / births_Total ; 
    
    tmp_EPrevalenceAmongPregnantWomen = ...
        (births_toHbEAgWomenLowVL + births_toHbEAgWomenHighVL) / (births_toHbSAgWomenLowVL + births_toHbSAgWomenHighVL + births_toHbEAgWomenLowVL + births_toHbEAgWomenHighVL) ;     
    
    % Mean year of birth of pregnant women
    tmp_MeanYearOfBirthOfPregnantWomen = time - ... 
        (sum(ages .* fert' .* squeeze(sum(sum(X([1:10 12 13:15],:,1,:),1),4))) / ...
            sum(fert' .* squeeze(sum(sum(X([1:10 12 13:15],:,1,:),1),4))));

      
    % increment the pointer
    itt = itt + 1;
    % increases every 0.1 years
    
end % end "time = TimeSteps" for loop



index_first_year_output=find(Time>=first_year_output, 1);
index_last_year_output=find(Time>=last_year_output, 1);

output.Time = Time(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.PrevEAg = PrevEAg_of_SAg_5yr(:,:,index_first_year_output:index_last_year_output); % 2 x 20 x num_years_output
output.NewChronicInfectionRate = Incid_chronic_all_5yr_approx(:,:,index_first_year_output:index_last_year_output); % 2 x 20 x num_years_output
output.Tot_Pop_1yr = Tot_Pop_1yr(:,:,index_first_year_output:index_last_year_output); % 2 x 100 x num_years_output
output.NewChronicInfectionRate_NeonatesOnly = Incid_babies_chronic_1yr_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.NumDecompCirr = NumDecompCirr(index_first_year_output:index_last_year_output); % 1 x num_years_output
assert(max(abs(...
    squeeze(sum(sum(Prev_Decomp_Cirr_1yr(:,:,index_first_year_output:index_last_year_output),1),2)) - ...
    NumDecompCirr(index_first_year_output:index_last_year_output)'...
    )) < 1e-9); % squeeze(sum(sum(Prev_Decomp_Cirr_1yr,1),2)) is a num_years_output x 1 matrix
output.NumLiverCancer = NumLiverCancer(index_first_year_output:index_last_year_output); % 1 x num_years_output
assert(max(abs(...
    squeeze(sum(sum(Prev_Liver_Cancer_1yr(:,:,index_first_year_output:index_last_year_output),1),2)) - ...
    NumLiverCancer(index_first_year_output:index_last_year_output)'...
    )) < 1e-8); % squeeze(sum(sum(Prev_Liver_Cancer_1yr,1),2)) is a num_years_output x 1 matrix
output.NumSAg_1yr = NumSAg_1yr(:,:,index_first_year_output:index_last_year_output); % 2 x 100 x num_years_output
output.NumEAg_chronic_1yr = NumEAg_chronic_1yr(:,:,index_first_year_output:index_last_year_output); % 2 x 100 x num_years_output
output.NumEAg_chronic_acute_1yr = NumEAg_chronic_acute_1yr(:,:,index_first_year_output:index_last_year_output); % 2 x 100 x num_years_output
output.yld_1yr = yld_1yr(:,:,index_first_year_output:index_last_year_output); % 2 x 100 x num_years_output
output.Incid_Deaths_1yr_approx = Incid_Deaths_1yr_approx(:,:,index_first_year_output:index_last_year_output); % 2 x 100 x num_years_output
output.Prev_Deaths_1yr = Prev_Deaths_1yr(:,:,index_first_year_output:index_last_year_output); % 2 x 100 x num_years_output
output.num_births_toHbEAgWomenHVL_1yr_approx = num_births_toHbEAgWomenHVL_1yr_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.num_births_toHbEAgWomenLVL_1yr_approx = num_births_toHbEAgWomenLVL_1yr_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.num_births_toHbSAgWomenHVL_1yr_approx = num_births_toHbSAgWomenHVL_1yr_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.num_births_toHbSAgWomenLVL_1yr_approx = num_births_toHbSAgWomenLVL_1yr_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.num_births_1yr_approx = num_births_1yr_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.num_births_chronic_HbEAgWomenHVL_1yr_approx = num_births_chronic_HbEAgWomenHVL_1yr_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.num_births_chronic_HbEAgWomenLVL_1yr_approx = num_births_chronic_HbEAgWomenLVL_1yr_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.num_births_chronic_HbSAgWomenHVL_1yr_approx = num_births_chronic_HbSAgWomenHVL_1yr_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.num_births_chronic_HbSAgWomenLVL_1yr_approx = num_births_chronic_HbSAgWomenLVL_1yr_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.RateBirthDoseVacc = RateBirthDoseVacc(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.RateInfantVacc = RateInfantVacc(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.PeripartumTreatment_HbEAg_HighVL_approx = PeripartumTreatment_HbEAg_HighVL_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.PeripartumTreatment_HbEAg_LowVL_approx = PeripartumTreatment_HbEAg_LowVL_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.PeripartumTreatment_HbSAg_HighVL_approx = PeripartumTreatment_HbSAg_HighVL_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.PeripartumTreatment_HbSAg_LowVL_approx = PeripartumTreatment_HbSAg_LowVL_approx(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.RatePeripartumTreatment = RatePeripartumTreatment(index_first_year_output:index_last_year_output); % 1 x num_years_output
output.PregnantWomenNeedToScreen = PregnantWomenNeedToScreen(index_first_year_output:index_last_year_output); % 1 x num_years_output; added 13.9.15
output.HBVPregnantWomenNeedToEvaluate = HBVPregnantWomenNeedToEvaluate(index_first_year_output:index_last_year_output); % 1 x num_years_output

output.beta_U5 = beta_U5;
output.p_HbSAg_av = p_HbSAg_av;

output.p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL = p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL;
output.p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL = p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL;

output.p_VerticalTransmission_HbSAgLowVL_NoIntv = p_VerticalTransmission_HbSAgLowVL_NoIntv;
output.p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc = p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc;
output.p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP = p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP;
output.p_VerticalTransmission_HbSAgLowVL_PAP = p_VerticalTransmission_HbSAgLowVL_PAP;
output.p_VerticalTransmission_HbSAgHighVL_NoIntv = p_VerticalTransmission_HbSAgHighVL_NoIntv;
output.p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc = p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc;
output.p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP = p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP;
output.p_VerticalTransmission_HbSAgHighVL_PAP = p_VerticalTransmission_HbSAgHighVL_PAP;

output.p_VerticalTransmission_HbEAgLowVL_NoIntv = p_VerticalTransmission_HbEAgLowVL_NoIntv;
output.p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc = p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc;
output.p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP = p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP;
output.p_VerticalTransmission_HbEAgLowVL_PAP = p_VerticalTransmission_HbEAgLowVL_PAP;
output.p_VerticalTransmission_HbEAgHighVL_NoIntv = p_VerticalTransmission_HbEAgHighVL_NoIntv;
output.p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc = p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc;
output.p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP = p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP;
output.p_VerticalTransmission_HbEAgHighVL_PAP = p_VerticalTransmission_HbEAgHighVL_PAP;



outputs_nums_cell_array = {...
    'beta_U5',...
    'p_HbSAg_av',...
    'p_VerticalTransmission_HbSAg_NoIntv_Ratio_HighVL_to_LowVL',...
    'p_VerticalTransmission_HbEAg_NoIntv_Ratio_HighVL_to_LowVL',...
    'p_VerticalTransmission_HbSAgLowVL_NoIntv',...
    'p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc',...
    'p_VerticalTransmission_HbSAgLowVL_BirthDoseVacc_PAP',...
    'p_VerticalTransmission_HbSAgLowVL_PAP',...
    'p_VerticalTransmission_HbSAgHighVL_NoIntv',...
    'p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc',...
    'p_VerticalTransmission_HbSAgHighVL_BirthDoseVacc_PAP',...
    'p_VerticalTransmission_HbSAgHighVL_PAP',...
    'p_VerticalTransmission_HbEAgLowVL_NoIntv',...
    'p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc',...
    'p_VerticalTransmission_HbEAgLowVL_BirthDoseVacc_PAP',...
    'p_VerticalTransmission_HbEAgLowVL_PAP',...
    'p_VerticalTransmission_HbEAgHighVL_NoIntv',...
    'p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc',...
    'p_VerticalTransmission_HbEAgHighVL_BirthDoseVacc_PAP',...
    'p_VerticalTransmission_HbEAgHighVL_PAP'...
    };
num_outputs_nums = length(outputs_nums_cell_array);
outputs_vectors_cell_array = {...
    'Time',...
    'NewChronicInfectionRate_NeonatesOnly',...
    'NumDecompCirr',...
    'NumLiverCancer',...
    'num_births_toHbEAgWomenHVL_1yr_approx',...
    'num_births_toHbEAgWomenLVL_1yr_approx',...
    'num_births_toHbSAgWomenHVL_1yr_approx',...
    'num_births_toHbSAgWomenLVL_1yr_approx',...
    'num_births_1yr_approx',...
    'num_births_chronic_HbEAgWomenHVL_1yr_approx',...
    'num_births_chronic_HbEAgWomenLVL_1yr_approx',...
    'num_births_chronic_HbSAgWomenHVL_1yr_approx',...
    'num_births_chronic_HbSAgWomenLVL_1yr_approx',...
    'RateBirthDoseVacc',...
    'RateInfantVacc',...
    'PeripartumTreatment_HbEAg_HighVL_approx',...
    'PeripartumTreatment_HbEAg_LowVL_approx',...
    'PeripartumTreatment_HbSAg_HighVL_approx',...
    'PeripartumTreatment_HbSAg_LowVL_approx',...
    'RatePeripartumTreatment',...
    'PregnantWomenNeedToScreen',...
    'HBVPregnantWomenNeedToEvaluate'...
    };
num_outputs_vectors = length(outputs_vectors_cell_array);
outputs_3D_cell_array = {...
    'PrevEAg',...
    'NewChronicInfectionRate',...
    'Tot_Pop_1yr',...
    'NumSAg_1yr',...
    'NumEAg_chronic_1yr',...
    'NumEAg_chronic_acute_1yr',...
    'yld_1yr',...
    'Incid_Deaths_1yr_approx',...
    'Prev_Deaths_1yr'...
    };
num_outputs_3D = length(outputs_3D_cell_array);
assert(all(ismember(fields(output),[outputs_nums_cell_array,outputs_vectors_cell_array,outputs_3D_cell_array])))
assert(all(ismember([outputs_nums_cell_array,outputs_vectors_cell_array,outputs_3D_cell_array],fields(output)))) % cell arrays contain the same elements
assert(isequal(sort(fields(output)),sort([outputs_nums_cell_array,outputs_vectors_cell_array,outputs_3D_cell_array]')))

for ii=1:num_outputs_nums
    fieldname = outputs_nums_cell_array{ii};
    field_output = output.(fieldname);
    assert(length(field_output)==1)
end

for ii=1:num_outputs_vectors
    fieldname = outputs_vectors_cell_array{ii};
    field_output = output.(fieldname);
    assert(length(field_output)==num_years_output)
end

for ii=1:num_outputs_3D
    fieldname = outputs_3D_cell_array{ii};
    field_output = output.(fieldname);
    assert(size(field_output,3)==num_years_output)
end





% Diagnositic work-up required for this simulation:
output.Dx_At_ANC_HBsAG = 0;
output.Dx_At_ANC_HBeAG = 0;
output.Dx_At_ANC_VL = 0; 

% Consider PAP for those who get BD:
if  (0==cov_BirthDoseAndTDF_EAgHighVL) && (0==cov_BirthDoseAndTDF_SAgHighVL) && (0==cov_BirthDoseAndTDF_EAgLowVL) && (0==cov_BirthDoseAndTDF_SAgLowVL) 
    
        % Zero coverage of all types of PAP 
        output.Dx_At_ANC_HBsAG = 0;
        output.Dx_At_ANC_HBeAG = 0;
        output.Dx_At_ANC_VL = 0; 
     
elseif (   (cov_BirthDoseAndTDF_EAgHighVL==cov_BirthDoseAndTDF_SAgHighVL) && (cov_BirthDoseAndTDF_SAgHighVL==cov_BirthDoseAndTDF_EAgLowVL) ...
                && (cov_BirthDoseAndTDF_EAgLowVL==cov_BirthDoseAndTDF_SAgLowVL) && (cov_BirthDoseAndTDF_SAgLowVL==cov_BirthDoseAndTDF_SAgHighVL)  )
        % PAP not differentiated by E or VL, so just use HBSAG
        output.Dx_At_ANC_HBsAG = 1;
        output.Dx_At_ANC_HBeAG = 0;
        output.Dx_At_ANC_VL = 0; 
    
elseif (   (cov_BirthDoseAndTDF_EAgHighVL==cov_BirthDoseAndTDF_EAgLowVL) && (cov_BirthDoseAndTDF_SAgHighVL==cov_BirthDoseAndTDF_SAgLowVL) ...
                && (cov_BirthDoseAndTDF_SAgLowVL~=cov_BirthDoseAndTDF_EAgLowVL) && (cov_BirthDoseAndTDF_SAgHighVL~=cov_BirthDoseAndTDF_EAgHighVL)    )    
        % PAP is differentiated only by E/S
        output.Dx_At_ANC_HBsAG = 1;
        output.Dx_At_ANC_HBeAG = 1;
        output.Dx_At_ANC_VL = 0;  
     
     
elseif (   (cov_BirthDoseAndTDF_SAgHighVL==cov_BirthDoseAndTDF_EAgHighVL) && (cov_BirthDoseAndTDF_SAgLowVL==cov_BirthDoseAndTDF_EAgLowVL) ...
                && (cov_BirthDoseAndTDF_SAgLowVL~=cov_BirthDoseAndTDF_SAgHighVL) && (cov_BirthDoseAndTDF_EAgLowVL~=cov_BirthDoseAndTDF_EAgHighVL) )
         % PAP is differentiated only by VL   
        output.Dx_At_ANC_HBsAG = 1;
        output.Dx_At_ANC_HBeAG = 0;
        output.Dx_At_ANC_VL = 1;         

elseif (   (cov_BirthDoseAndTDF_EAgHighVL~=cov_BirthDoseAndTDF_SAgHighVL) && (cov_BirthDoseAndTDF_SAgHighVL~=cov_BirthDoseAndTDF_EAgLowVL) ...
                && (cov_BirthDoseAndTDF_EAgLowVL~=cov_BirthDoseAndTDF_SAgLowVL) && (cov_BirthDoseAndTDF_SAgLowVL~=cov_BirthDoseAndTDF_SAgHighVL)  )
        % PAP is differentiated by both E and VL
        output.Dx_At_ANC_HBsAG = 1;
        output.Dx_At_ANC_HBeAG = 1;
        output.Dx_At_ANC_VL = 1;    
         
else
        % Fail
        assert(false)
         
end
    

% Consider those who do not get BD
if (0==cov_TDFOnly_EAgHighVL) && (0==cov_TDFOnly_SAgHighVL) && (0==cov_TDFOnly_EAgLowVL) && (0==cov_TDFOnly_SAgLowVL) 
        % Zero coverage so do nothing
    
elseif (   (cov_TDFOnly_EAgHighVL==cov_TDFOnly_SAgHighVL) && (cov_TDFOnly_SAgHighVL==cov_TDFOnly_EAgLowVL) ...
                && (cov_TDFOnly_EAgLowVL==cov_TDFOnly_SAgLowVL) && (cov_TDFOnly_SAgLowVL==cov_TDFOnly_EAgHighVL) )
        % PAP not differentiated by E or VL, so just use HBSAG
        output.Dx_At_ANC_HBsAG = min(1,output.Dx_At_ANC_HBsAG+1);
    
elseif (   (cov_TDFOnly_EAgHighVL==cov_TDFOnly_EAgLowVL) && (cov_TDFOnly_SAgHighVL==cov_TDFOnly_SAgLowVL) ...
                && (cov_TDFOnly_SAgLowVL~=cov_TDFOnly_EAgLowVL) && (cov_TDFOnly_SAgHighVL~=cov_TDFOnly_EAgHighVL)     )
        % PAP is differentiated only by E/S
        output.Dx_At_ANC_HBsAG = min(1,output.Dx_At_ANC_HBsAG+1);
        output.Dx_At_ANC_HBeAG = min(1,output.Dx_At_ANC_HBeAG+1);
     
elseif (   (cov_TDFOnly_SAgHighVL==cov_TDFOnly_SAgHighVL) && (cov_TDFOnly_SAgLowVL==cov_TDFOnly_SAgLowVL) ...
                && (cov_TDFOnly_SAgLowVL~=cov_TDFOnly_SAgHighVL) && (cov_TDFOnly_EAgLowVL~=cov_TDFOnly_EAgHighVL)  )
        % PAP is differentiated only by VL   
        output.Dx_At_ANC_HBsAG = min(1,output.Dx_At_ANC_HBsAG+1);
        output.Dx_At_ANC_VL = min(1,output.Dx_At_ANC_VL+1);        

elseif (   (cov_TDFOnly_EAgHighVL==cov_TDFOnly_SAgHighVL) && (cov_TDFOnly_SAgHighVL==cov_TDFOnly_EAgLowVL) ...
                && (cov_TDFOnly_EAgLowVL==cov_TDFOnly_SAgLowVL) && (cov_TDFOnly_SAgLowVL==cov_TDFOnly_SAgHighVL) )
        % PAP is differentiated by both E and VL
        output.Dx_At_ANC_HBsAG = min(1,output.Dx_At_ANC_HBsAG+1);
        output.Dx_At_ANC_HBeAG = min(1,output.Dx_At_ANC_HBeAG+1);
        output.Dx_At_ANC_VL = min(1,output.Dx_At_ANC_VL+1);    
    
end




assert(length(output.Dx_At_ANC_HBsAG)==1)
assert(length(output.Dx_At_ANC_HBeAG)==1)
assert(length(output.Dx_At_ANC_VL)==1)





fields_of_interest = fields(output);
num_fields = length(fields_of_interest);

for ii=1:num_fields
    curr_field = fields_of_interest{ii};
    assert(all(all(all(output.(curr_field)>=0))))
end



end % end function HBVmodel





function cov_vec = make_coverage_vector(xvals_vec,yvals_vec,scaleup_years_vec,TimeSteps)
% TimeSteps contains all of the time steps in the run (from 1890 to the last year of the run)
% yvals_vec contains important coverage values from 1980 to the last year of the run
% xvals_vec contains the years that correspond to the values in yvals_vec
% scaleup_years_vec contains the years between which linear scale-up in coverage occurs

    % checks
    assert(isequal(size(xvals_vec),size(yvals_vec)))
    assert(isequal(xvals_vec,sort(xvals_vec))) % ensure that years are in order
    assert(isequal(xvals_vec,unique(xvals_vec))) % ensure that there are no duplicate years
    assert(xvals_vec(1)>=TimeSteps(1))
    assert(xvals_vec(end)<=TimeSteps(end))
    assert(yvals_vec(1)==0)
    assert(isequal(size(scaleup_years_vec),[1 2])) % first and last year of scale-up
    assert(scaleup_years_vec(1)<scaleup_years_vec(end))
    assert(all(ismember(scaleup_years_vec,xvals_vec)))

    % set up vectors
    years_vec = TimeSteps;
    [~,cov_vec_indices] = ismember(xvals_vec,years_vec); % [Lia,Locb] = ismember(A,B): Locb contains the lowest index in B for each value in A that is a member of B
    assert(all(cov_vec_indices)>0) % Values of 0 indicate where A is not a member of B.
    assert(isequal(size(xvals_vec),size(cov_vec_indices)))

    % assign values to cov_vec
    cov_vec = -99*ones(1,length(years_vec));
    last_pos = cov_vec_indices(1);
    num_vals_needed = length(years_vec(1:last_pos)); % number of coverage values needed between the first year in years_vec and first year in xvals_vec
    cov_vec(1:cov_vec_indices(1)) = repmat(yvals_vec(1),1,num_vals_needed);
    num_gaps = length(xvals_vec) - 1;
    for gap_num=1:num_gaps
        first_pos = cov_vec_indices(gap_num); % position of lefthand number in years_vec
        last_pos = cov_vec_indices(gap_num+1); % position of righthand number in years_vec
        assert(first_pos<last_pos)
        num_vals_needed = length(years_vec(first_pos:last_pos)) - 1; % number of coverage values needed between successive years in xvals_vec
        cov_vec(first_pos:last_pos) = [repmat(yvals_vec(gap_num),1,num_vals_needed) yvals_vec(gap_num+1)];
    end
    num_vals_needed = length(years_vec(last_pos:end)); % number of coverage values needed between last year in xvals_vec and the last year in years_vec
    cov_vec(last_pos:end) = repmat(yvals_vec(gap_num+1),1,num_vals_needed);

    % assign scale-up values to cov_vec
    first_scaleup_pos_in_yvals_vec = find(scaleup_years_vec(1)==xvals_vec);
    last_scaleup_pos_in_yvals_vec = find(scaleup_years_vec(end)==xvals_vec);
    assert(first_scaleup_pos_in_yvals_vec+1==last_scaleup_pos_in_yvals_vec)
    first_scaleup_pos_in_years_vec = find(scaleup_years_vec(1)==years_vec);
    last_scaleup_pos_in_years_vec = find(scaleup_years_vec(end)==years_vec);
    assert(first_scaleup_pos_in_years_vec<last_scaleup_pos_in_years_vec)
    cov_vec(first_scaleup_pos_in_years_vec:last_scaleup_pos_in_years_vec) = interp1(...
        [scaleup_years_vec(1) scaleup_years_vec(end)],...
        [yvals_vec(first_scaleup_pos_in_yvals_vec) yvals_vec(last_scaleup_pos_in_yvals_vec)],...
        years_vec(first_scaleup_pos_in_years_vec:last_scaleup_pos_in_years_vec),...
        'linear');

    % checks
    assert(all(cov_vec)>=0) % ensure that there are no -99 values left
    assert(all(cov_vec)<=1)
    assert(isequal(size(cov_vec),size(years_vec))) % ensure that cov_vec has not been extended in length

end % end make_coverage_vector function

