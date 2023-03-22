function Costs_map_countries = cost_model_countries(Runs_map, screen_every_pregnancy, do_ESLD, ...
    UnitCost_struct, AverageBabiesPerWoman_map, num_years_1980_2100, num_strategies, strategies_vec, num_stochas_runs, num_parameter_values)


    if exist('num_parameter_values')~=1 % i.e. is not a variable in the workspace
        num_parameter_values = 1;
    end

    if num_stochas_runs>1
        assert(num_parameter_values==1)
    end


    assert((UnitCost_struct.UnitCost_ScreeningANC_LOW<UnitCost_struct.UnitCost_ScreeningANC) && (UnitCost_struct.UnitCost_ScreeningANC<UnitCost_struct.UnitCost_ScreeningANC_HIGH))
    assert((UnitCost_struct.UnitCost_HBeAg_LOW<UnitCost_struct.UnitCost_HBeAg) && (UnitCost_struct.UnitCost_HBeAg<UnitCost_struct.UnitCost_HBeAg_HIGH))
    assert((UnitCost_struct.UnitCost_VL_LOW<UnitCost_struct.UnitCost_VL) && (UnitCost_struct.UnitCost_VL<UnitCost_struct.UnitCost_VL_HIGH))
    assert((UnitCost_struct.UnitCost_PPT_Drug_Low==UnitCost_struct.UnitCost_PPT_Drug) && (UnitCost_struct.UnitCost_PPT_Drug==UnitCost_struct.UnitCost_PPT_Drug_High))
    assert((UnitCost_struct.UnitCost_PPT_Monitoring_Low<UnitCost_struct.UnitCost_PPT_Monitoring) && (UnitCost_struct.UnitCost_PPT_Monitoring<UnitCost_struct.UnitCost_PPT_Monitoring_High))
    assert((UnitCost_struct.UnitCost_HBIG_Low==0) && (UnitCost_struct.UnitCost_HBIG==0) && (UnitCost_struct.UnitCost_HBIG_High==0))
    if do_ESLD
        assert(UnitCost_struct.UnitCost_DC>0)
        assert(UnitCost_struct.UnitCost_DC_Low>0)
        assert(UnitCost_struct.UnitCost_DC_High==0)
        assert(UnitCost_struct.UnitCost_HCC>0)
        assert(UnitCost_struct.UnitCost_HCC_Low>0)
        assert(UnitCost_struct.UnitCost_HCC_High==0)
    else
        assert(UnitCost_struct.UnitCost_DC==0)
        assert(UnitCost_struct.UnitCost_DC_Low==0)
        assert(UnitCost_struct.UnitCost_DC_High==0)
        assert(UnitCost_struct.UnitCost_HCC==0)
        assert(UnitCost_struct.UnitCost_HCC_Low==0)
        assert(UnitCost_struct.UnitCost_HCC_High==0)
    end


    if num_parameter_values>1
        assert(num_stochas_runs==1)
        expand_matrix_func = @(xxx) repmat(reshape(xxx,1,1,[]),num_parameter_values,num_years_1980_2100); % expand a vector into a 3D matrix
    else
        expand_matrix_func = @(xxx) repmat(reshape(xxx,1,1,[]),num_stochas_runs,num_years_1980_2100); % expand a vector into a 3D matrix
    end
    UnitCost_HBIG_mat = expand_matrix_func([UnitCost_struct.UnitCost_HBIG_Low UnitCost_struct.UnitCost_HBIG UnitCost_struct.UnitCost_HBIG_High]);
    UnitCost_PPT_mat = expand_matrix_func(...
        [UnitCost_struct.UnitCost_PPT_Drug_Low UnitCost_struct.UnitCost_PPT_Drug UnitCost_struct.UnitCost_PPT_Drug_High] + ...
        [UnitCost_struct.UnitCost_PPT_Monitoring_Low UnitCost_struct.UnitCost_PPT_Monitoring UnitCost_struct.UnitCost_PPT_Monitoring_High]...
        );
    UnitCost_ANCScreening_mat = expand_matrix_func([UnitCost_struct.UnitCost_ScreeningANC_LOW UnitCost_struct.UnitCost_ScreeningANC UnitCost_struct.UnitCost_ScreeningANC_HIGH]);
    UnitCost_HBeAg_mat = expand_matrix_func([UnitCost_struct.UnitCost_HBeAg_LOW UnitCost_struct.UnitCost_HBeAg UnitCost_struct.UnitCost_HBeAg_HIGH]);
    UnitCost_VL_mat = expand_matrix_func([UnitCost_struct.UnitCost_VL_LOW UnitCost_struct.UnitCost_VL UnitCost_struct.UnitCost_VL_HIGH]);
    UnitCost_DC_mat = expand_matrix_func([UnitCost_struct.UnitCost_DC_Low UnitCost_struct.UnitCost_DC UnitCost_struct.UnitCost_DC_High]);
    UnitCost_HCC_mat = expand_matrix_func([UnitCost_struct.UnitCost_HCC_Low UnitCost_struct.UnitCost_HCC UnitCost_struct.UnitCost_HCC_High]);
    if num_parameter_values>1
        assert(isequal(size(UnitCost_HBIG_mat),[num_parameter_values num_years_1980_2100 3]))
        assert(isequal(size(UnitCost_PPT_mat),[num_parameter_values num_years_1980_2100 3]))
        assert(isequal(size(UnitCost_ANCScreening_mat),[num_parameter_values num_years_1980_2100 3]))
        assert(isequal(size(UnitCost_HBeAg_mat),[num_parameter_values num_years_1980_2100 3]))
        assert(isequal(size(UnitCost_VL_mat),[num_parameter_values num_years_1980_2100 3]))
        assert(isequal(size(UnitCost_DC_mat),[num_parameter_values num_years_1980_2100 3]))
        assert(isequal(size(UnitCost_HCC_mat),[num_parameter_values num_years_1980_2100 3]))
    else
        assert(isequal(size(UnitCost_HBIG_mat),[num_stochas_runs num_years_1980_2100 3]))
        assert(isequal(size(UnitCost_PPT_mat),[num_stochas_runs num_years_1980_2100 3]))
        assert(isequal(size(UnitCost_ANCScreening_mat),[num_stochas_runs num_years_1980_2100 3]))
        assert(isequal(size(UnitCost_HBeAg_mat),[num_stochas_runs num_years_1980_2100 3]))
        assert(isequal(size(UnitCost_VL_mat),[num_stochas_runs num_years_1980_2100 3]))
        assert(isequal(size(UnitCost_DC_mat),[num_stochas_runs num_years_1980_2100 3]))
        assert(isequal(size(UnitCost_HCC_mat),[num_stochas_runs num_years_1980_2100 3]))
    end


    country_list = keys(Runs_map);
    num_countries = length(country_list);


    Costs_map_countries = containers.Map;


    for country_num=1:num_countries


        ISO = country_list{country_num};
        countries_runs_stochas_cell_array = Runs_map(ISO);
        assert(size(countries_runs_stochas_cell_array,2)==num_strategies)


        AverageBabiesPerWoman = AverageBabiesPerWoman_map(ISO);
        AverageBabiesPerWomanLow = AverageBabiesPerWoman/2;
        if screen_every_pregnancy
            ANCScreening_mat = UnitCost_ANCScreening_mat;
        else
            if do_ESLD
                AverageBabiesPerWoman_mat = expand_matrix_func(repmat(AverageBabiesPerWoman,1,3));
            else
                AverageBabiesPerWoman_mat = expand_matrix_func([AverageBabiesPerWoman AverageBabiesPerWoman AverageBabiesPerWomanLow]);
            end
            assert(isequal(size(UnitCost_ANCScreening_mat),size(AverageBabiesPerWoman_mat)))
            ANCScreening_mat = UnitCost_ANCScreening_mat ./ AverageBabiesPerWoman_mat;
        end


        Runs_1 = countries_runs_stochas_cell_array{1};
        if num_parameter_values>1
            assert(isequal(size(Runs_1.RateInfantVacc),[num_parameter_values num_years_1980_2100]))
            assert(isequal(size(Runs_1.RateBirthDoseVacc),[num_parameter_values num_years_1980_2100]))
            assert(isequal(size(Runs_1.PregnantWomenNeedToScreen),[num_parameter_values num_years_1980_2100]))
            assert(isequal(size(Runs_1.HBVPregnantWomenNeedToEvaluate),[num_parameter_values num_years_1980_2100]))
            assert(isequal(size(Runs_1.RatePeripartumTreatment),[num_parameter_values num_years_1980_2100]))
            assert(isequal(size(Runs_1.NumDecompCirr),[num_parameter_values num_years_1980_2100]))
            assert(isequal(size(Runs_1.NumLiverCancer),[num_parameter_values num_years_1980_2100]))
        else
            assert(isequal(size(Runs_1.RateInfantVacc),[num_stochas_runs num_years_1980_2100]))
            assert(isequal(size(Runs_1.RateBirthDoseVacc),[num_stochas_runs num_years_1980_2100]))
            assert(isequal(size(Runs_1.PregnantWomenNeedToScreen),[num_stochas_runs num_years_1980_2100]))
            assert(isequal(size(Runs_1.HBVPregnantWomenNeedToEvaluate),[num_stochas_runs num_years_1980_2100]))
            assert(isequal(size(Runs_1.RatePeripartumTreatment),[num_stochas_runs num_years_1980_2100]))
        end
        assert(isequal(size(Runs_1.Dx_At_ANC_HBsAG),[1 1]))
        assert(isequal(size(Runs_1.Dx_At_ANC_HBeAG),[1 1]))
        assert(isequal(size(Runs_1.Dx_At_ANC_VL),[1 1]))


        %% APPLY COSTING MODEL TO THE MODEL RUN APPLIED
        % Each Cost Category is given as a 3-column set.
        %   First column: low-bound estimate
        %   Second column: best estimate
        %   Third column: high-bound estimate

        country_costs_cell_array = cell(1,num_strategies);

        for strategy_num=strategies_vec

            output = [];

            output.Cost_InfantVacc = repmat(UnitCost_struct.UnitCost_InfantVacc * countries_runs_stochas_cell_array{strategy_num}.RateInfantVacc,[1 1 3]);

            output.Cost_BD = repmat(UnitCost_struct.UnitCost_BirthDoseVacc * countries_runs_stochas_cell_array{strategy_num}.RateBirthDoseVacc,[1 1 3]);

            output.Cost_HBIG = zeros(size(output.Cost_InfantVacc));
            if any(UnitCost_HBIG_mat)>0
                output.Cost_HBIG = UnitCost_HBIG_mat .* repmat(countries_runs_stochas_cell_array{strategy_num}.HBVPregnantWomenNeedToEvaluate,[1 1 3]);
            end

            % Cost of PPT: including drug and monitoring costs 
            output.Cost_PPT = UnitCost_PPT_mat .* repmat(countries_runs_stochas_cell_array{strategy_num}.RatePeripartumTreatment,[1 1 3]);

            % Cost of the diagnostics:
            output.Cost_ANCScreening = zeros(size(output.Cost_InfantVacc));
            if countries_runs_stochas_cell_array{strategy_num}.Dx_At_ANC_HBsAG>0     % test to see if any screening is done for pregnant women
                output.Cost_ANCScreening = ANCScreening_mat .* repmat(countries_runs_stochas_cell_array{strategy_num}.PregnantWomenNeedToScreen,[1 1 3]);
            end

            % Cost of HBeAg test for all preg women who test HBsAg positive
            output.Cost_HBeAg = zeros(size(output.Cost_InfantVacc));
            if countries_runs_stochas_cell_array{strategy_num}.Dx_At_ANC_HBeAG>0     % Test to see if e-ag screening is used (if not, zero costs for e-ag)
                output.Cost_HBeAg = UnitCost_HBeAg_mat .* repmat(countries_runs_stochas_cell_array{strategy_num}.HBVPregnantWomenNeedToEvaluate,[1 1 3]);
            end

            % Cost of VL test for all preg women tested HBsAg positive
            output.Cost_VL = zeros(size(output.Cost_InfantVacc));
            if countries_runs_stochas_cell_array{strategy_num}.Dx_At_ANC_VL>0
                output.Cost_VL = UnitCost_VL_mat .* repmat(countries_runs_stochas_cell_array{strategy_num}.HBVPregnantWomenNeedToEvaluate,[1 1 3]);
            end

            % Costs of management of ESLD
            % DC
            output.Cost_DC = zeros(size(output.Cost_InfantVacc));
            if any(UnitCost_DC_mat)>0
                output.Cost_DC = UnitCost_DC_mat .* repmat(countries_runs_stochas_cell_array{strategy_num}.NumDecompCirr,[1 1 3]);
            end

            % Costs of management of ESLD
            % HCC                                         
            output.Cost_HCC = zeros(size(output.Cost_InfantVacc));
            if any(UnitCost_HCC_mat)>0
                output.Cost_HCC = UnitCost_HCC_mat .* repmat(countries_runs_stochas_cell_array{strategy_num}.NumLiverCancer,[1 1 3]);
            end

            output.Total = ...
                output.Cost_InfantVacc + ...
                output.Cost_BD + ...
                output.Cost_PPT + ...
                output.Cost_ANCScreening + ...
                output.Cost_HBIG + ...
                output.Cost_HBeAg + ...
                output.Cost_VL + ...
                output.Cost_HCC + ...
                output.Cost_DC;


            country_costs_cell_array{strategy_num} = output;
            
            
        end % end strategy_num for loop

        clear strategy_num
        assert(length(country_costs_cell_array)==num_strategies)
        if num_parameter_values>1
            assert(isequal(size(country_costs_cell_array{1}.Cost_InfantVacc),[num_parameter_values num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_BD),[num_parameter_values num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_HBIG),[num_parameter_values num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_PPT),[num_parameter_values num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_ANCScreening),[num_parameter_values num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_HBeAg),[num_parameter_values num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_VL),[num_parameter_values num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_DC),[num_parameter_values num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_HCC),[num_parameter_values num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Total),[num_parameter_values num_years_1980_2100 3]))
        else
            assert(isequal(size(country_costs_cell_array{1}.Cost_InfantVacc),[num_stochas_runs num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_BD),[num_stochas_runs num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_HBIG),[num_stochas_runs num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_PPT),[num_stochas_runs num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_ANCScreening),[num_stochas_runs num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_HBeAg),[num_stochas_runs num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_VL),[num_stochas_runs num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_DC),[num_stochas_runs num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Cost_HCC),[num_stochas_runs num_years_1980_2100 3]))
            assert(isequal(size(country_costs_cell_array{1}.Total),[num_stochas_runs num_years_1980_2100 3]))
        end


        Costs_map_countries(ISO) = country_costs_cell_array;

    
    end % end country_num for loop


    assert(length(Costs_map_countries)==num_countries)


end % end function cost_model_countries


