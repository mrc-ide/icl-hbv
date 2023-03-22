function out = get_range_of_FracSPosHighVL(FracEPosHighVL,p_HbSAg_av, p_HbEAg_av)

% For a given value of:
%   * - FracEPosHighVL (Fraction of e-positive that are high VL)
%   * - p_HbSAg_av (Average transmission for s-positive, aX VL states)
%   * - p_HbEAg_av (Average transmission for e-positive, aX VL states)
% ... find the high and low values for FracSPosHighVL that satisfy the
% constraint of implied transmission rates (by e/s and VL) being the same
% irrespective of E/S status and being bounded [0,1].


FracSPosHighVL_rng = linspace(0,1,1000);
AllowFlag = false(1,1000);

for i=1:length(FracSPosHighVL_rng)

    FracSPosHighVL = FracSPosHighVL_rng(i);
    
    beta_low_vl = max(1e-10,( p_HbEAg_av * FracSPosHighVL - p_HbSAg_av * FracEPosHighVL ) / ( FracSPosHighVL - FracEPosHighVL ) );
    beta_high_vl = (p_HbSAg_av - (1-FracSPosHighVL)*beta_low_vl)/FracSPosHighVL;

    
    if (beta_low_vl>0) && (beta_low_vl<1) && (beta_high_vl>0) && (beta_high_vl<1) && (beta_high_vl>beta_low_vl)
        AllowFlag(i) = true;
    end

end

allow_min = FracSPosHighVL_rng(find(AllowFlag,1,'first')); % the value in FracSPosHighVL_rng corresponding to the first index in AllowFlag containing a nonzero element
allow_max = FracSPosHighVL_rng(find(AllowFlag,1,'last')); % the value in FracSPosHighVL_rng corresponding to the last index in AllowFlag containing a nonzero element

% stem(FracSPosHighVL_rng,AllowFlag)
out = [allow_min, allow_max];

end


% 2nd equation: 
% beta_high_vl = (p_HbSAg_av - (1-FracSPosHighVL)*beta_low_vl)/FracSPosHighVL;
% beta_high_vl = p_HbSAg_av/FracSPosHighVL - ((1-FracSPosHighVL)*beta_low_vl)/FracSPosHighVL;
% p_HbSAg_av/FracSPosHighVL = beta_high_vl + ((1-FracSPosHighVL)*beta_low_vl)/FracSPosHighVL;
% p_HbSAg_av = FracSPosHighVL*beta_high_vl + (1-FracSPosHighVL)*beta_low_vl
% p_HbSAg_av = FracSPosHighVL*beta_high_vl + FracSPosLowVL*beta_low_vl  ***
% hence p_HbSAg_av is the weighted average of beta_high_vl and beta_low_vl,
% with the weights being the fraction of s+ people with high viral load and
% the fraction of s+ people with low viral load

% E equivalent of the 2nd equation:
% beta_high_vl = (p_HbEAg_av - (1-FracEPosHighVL)*beta_low_vl)/FracEPosHighVL;
% p_HbEAg_av = FracEPosHighVL*beta_high_vl + (1-FracEPosHighVL)*beta_low_vl  ***

% Hence, equating the S version with the E version gives:
% (p_HbEAg_av - (1-FracEPosHighVL)*beta_low_vl)/FracEPosHighVL = (p_HbSAg_av - (1-FracSPosHighVL)*beta_low_vl)/FracSPosHighVL
% 0 = (p_HbSAg_av - (1-FracSPosHighVL)*beta_low_vl)/FracSPosHighVL - (p_HbEAg_av - (1-FracEPosHighVL)*beta_low_vl)/FracEPosHighVL
% 0 = p_HbSAg_av/FracSPosHighVL - (1-FracSPosHighVL)*beta_low_vl/FracSPosHighVL - p_HbEAg_av/FracEPosHighVL + (1-FracEPosHighVL)*beta_low_vl/FracEPosHighVL
% (1-FracSPosHighVL)*beta_low_vl/FracSPosHighVL - (1-FracEPosHighVL)*beta_low_vl/FracEPosHighVL = p_HbSAg_av/FracSPosHighVL - p_HbEAg_av/FracEPosHighVL
% beta_low_vl*[(1-FracSPosHighVL)/FracSPosHighVL - (1-FracEPosHighVL)/FracEPosHighVL] = p_HbSAg_av/FracSPosHighVL - p_HbEAg_av/FracEPosHighVL
% beta_low_vl*[(1-FracSPosHighVL)*FracEPosHighVL/(FracSPosHighVL*FracEPosHighVL) - (1-FracEPosHighVL)*FracSPosHighVL/(FracEPosHighVL*FracSPosHighVL)] = p_HbSAg_av*FracEPosHighVL/(FracSPosHighVL*FracEPosHighVL) - p_HbEAg_av*FracSPosHighVL/(FracEPosHighVL*FracSPosHighVL)
% beta_low_vl*[((1-FracSPosHighVL)*FracEPosHighVL - (1-FracEPosHighVL)*FracSPosHighVL)/(FracEPosHighVL*FracSPosHighVL)] = (p_HbSAg_av*FracEPosHighVL - p_HbEAg_av*FracSPosHighVL)/(FracEPosHighVL*FracSPosHighVL)
% beta_low_vl = (p_HbSAg_av*FracEPosHighVL - p_HbEAg_av*FracSPosHighVL)/((1-FracSPosHighVL)*FracEPosHighVL - (1-FracEPosHighVL)*FracSPosHighVL)
% beta_low_vl = (p_HbEAg_av*FracSPosHighVL - p_HbSAg_av*FracEPosHighVL)/((1-FracEPosHighVL)*FracSPosHighVL - (1-FracSPosHighVL)*FracEPosHighVL)
% beta_low_vl = (p_HbEAg_av*FracSPosHighVL - p_HbSAg_av*FracEPosHighVL)/(FracSPosHighVL - FracEPosHighVL*FracSPosHighVL - FracEPosHighVL + FracSPosHighVL*FracEPosHighVL)
% beta_low_vl = (p_HbEAg_av*FracSPosHighVL - p_HbSAg_av*FracEPosHighVL)/(FracSPosHighVL - FracEPosHighVL)
% which is the 1st equation.

% Hence, all FracSPosHighVL_rng values for which beta_low_vl and 
% beta_high_vl are both between 0 and 1 and beta_low_vl < beta_high_vl, 
% where beta_high_vl and beta_low_vl are defined according the equations
% p_HbSAg_av = FracSPosHighVL*beta_high_vl + (1-FracSPosHighVL)*beta_low_vl
% and 
% p_HbEAg_av = FracEPosHighVL*beta_high_vl + (1-FracEPosHighVL)*beta_low_vl,
% are accepted.



