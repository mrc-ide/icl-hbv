function [UnitCost_struct, AverageBabiesPerWoman_map] = make_cost_objects_countries(do_ESLD,region_code_IHME_map,ISO_array)

    %% COST MODEL - ADAPTED FROM GLOBAL ANALYSIS
    % FROM SN labelled (updated 16.8.19)
    % added costs of ESLD (updated 12.11.19)




    % LOAD COST ASSUMPTIONS
    % All regions the same

    UnitCost_struct.UnitCost_InfantVacc = 1;
    UnitCost_struct.UnitCost_BirthDoseVacc = 1;

    UnitCost_struct.UnitCost_ScreeningANC = 1.6;            % HBsAg POC test cost
    UnitCost_struct.UnitCost_ScreeningANC_LOW = 0.4;
    UnitCost_struct.UnitCost_ScreeningANC_HIGH = 2.8;

    UnitCost_struct.UnitCost_HBeAg = 7.5;
    UnitCost_struct.UnitCost_HBeAg_LOW = 3;
    UnitCost_struct.UnitCost_HBeAg_HIGH = 40;

    UnitCost_struct.UnitCost_VL = 15;
    UnitCost_struct.UnitCost_VL_LOW = 5;
    UnitCost_struct.UnitCost_VL_HIGH = 100;

    UnitCost_struct.UnitCost_PPT_Drug = 10;
    UnitCost_struct.UnitCost_PPT_Drug_Low= 10;
    UnitCost_struct.UnitCost_PPT_Drug_High = 10;

    UnitCost_struct.UnitCost_PPT_Monitoring = 10;
    UnitCost_struct.UnitCost_PPT_Monitoring_Low = 5;
    UnitCost_struct.UnitCost_PPT_Monitoring_High = 40;

    UnitCost_struct.UnitCost_HBIG = 0;
    UnitCost_struct.UnitCost_HBIG_Low = 0;
    UnitCost_struct.UnitCost_HBIG_High = 0;

    if do_ESLD
        UnitCost_struct.UnitCost_ScreeningANC_LOW = UnitCost_struct.UnitCost_ScreeningANC;
        UnitCost_struct.UnitCost_ScreeningANC_HIGH = UnitCost_struct.UnitCost_ScreeningANC;

        UnitCost_struct.UnitCost_HBeAg_LOW = UnitCost_struct.UnitCost_HBeAg;
        UnitCost_struct.UnitCost_HBeAg_HIGH = UnitCost_struct.UnitCost_HBeAg;

        UnitCost_struct.UnitCost_VL_LOW = UnitCost_struct.UnitCost_VL;
        UnitCost_struct.UnitCost_VL_HIGH = UnitCost_struct.UnitCost_VL;

        UnitCost_struct.UnitCost_PPT_Monitoring_Low = UnitCost_struct.UnitCost_PPT_Monitoring;
        UnitCost_struct.UnitCost_PPT_Monitoring_High = UnitCost_struct.UnitCost_PPT_Monitoring;

        UnitCost_struct.UnitCost_DC = 500;
        UnitCost_struct.UnitCost_DC_Low = 2500;    %ie Highest costs averted -- > ie lowest ICER
        UnitCost_struct.UnitCost_DC_High = 0;      % BASECASE ie No costs averted -- > highest ICER

        UnitCost_struct.UnitCost_HCC = 500;
        UnitCost_struct.UnitCost_HCC_Low = 2500;   %ie Highest costs averted -- > ie lowest ICER
        UnitCost_struct.UnitCost_HCC_High = 0;    % BASECASE ie No costs averted -- > highest ICER
    else
        % For basecase, default value of ESLD = zero
        UnitCost_struct.UnitCost_DC = 0;
        UnitCost_struct.UnitCost_DC_Low = 0;    % ie Highest costs averted -- > ie lowest ICER
        UnitCost_struct.UnitCost_DC_High = 0;      % BASECASE ie No costs averted -- > highest ICER

        UnitCost_struct.UnitCost_HCC = 0;
        UnitCost_struct.UnitCost_HCC_Low = 0;   % ie Highest costs averted -- > ie lowest ICER
        UnitCost_struct.UnitCost_HCC_High = 0;     % BASECASE ie No costs averted -- > highest ICER
    end




    % Average Number of babies per woman (2010-2015) by region

    AverageBabiesPerWoman_vec = -99 * ones(1,21);

    % Asia
    AverageBabiesPerWoman_vec([1 2 3 5 18]) = 2.2;

    % Oceania
    AverageBabiesPerWoman_vec([4 19]) = 2.4;

    % Europe
    AverageBabiesPerWoman_vec([6 7 20]) = 1.6;

    % Africa
    AverageBabiesPerWoman_vec([8 9 10 11 12]) = 4.7;

    % LA & Carribean
    AverageBabiesPerWoman_vec([13 14 15 16 17]) = 2.2;

    % North America
    AverageBabiesPerWoman_vec(21) = 1.9;

    assert(all(AverageBabiesPerWoman_vec>0))

    AverageBabiesPerWoman_vals = cellfun(@(xx) AverageBabiesPerWoman_vec(region_code_IHME_map(xx)),ISO_array,'UniformOutput',false);
    AverageBabiesPerWoman_map = containers.Map(ISO_array,AverageBabiesPerWoman_vals);




end % end function make_cost_objects_countries

