%% SearchSaturationFactor

function [factor_k_new,...
          sol_full_new] = SearchSaturationFactor(model,mu,f,osenseStr,rxnID,...
                                                 f_transporter,kcat_glc,...
                                                 Info_enzyme,...
                                                 Info_mRNA,...
                                                 Info_protein,...
                                                 Info_ribosome,...
                                                 Info_tRNA,...
                                                 Check)
% do not consider glucose transporter, then choose 1
% consider glucose transporter, then choose 2
if exist('Check', 'var')
    if isempty(Check)
        Check = 1;
    end
else
    Check = 2;
end


if Check == 2

    factor_org = 1;

    fileName = WriteLPSatFactor(model,mu,f,osenseStr,rxnID,factor_org,...
                                f_transporter,kcat_glc,1,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
    command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
    system(command,'-echo');
    fileName_out = 'Simulation.lp.out';
    [~,sol_status,sol_full] = ReadSoplexResult(fileName_out,model);

    if strcmp(sol_status,'optimal')
        [~,~,f_enzyme_inact] = CheckInactiveEnzyme(model,sol_full);
    else 
        error('No solution is found for the original simulation.');
    end

    tmp_protein = f_enzyme_inact - f_transporter * 0.46;

    if tmp_protein <= 0
        factor_k_new = factor_org;
        sol_full_new = sol_full;
    else
        factor_high = 1;
        factor_low = 0;
    %     while factor_high-factor_low > 0.001 || tmp_protein > 0 || ~strcmp(sol_status,'optimal')
        while factor_high-factor_low > 0.0001 || (tmp_protein > 0 && factor_high-factor_low > 0.000001) || ~strcmp(sol_status,'optimal')
            factor_mid = (factor_low+factor_high)/2;
            disp(['kcat factor = ' num2str(factor_mid)]);

            fileName = WriteLPSatFactor(model,mu,f,osenseStr,rxnID,factor_mid,...
                                        f_transporter,kcat_glc*10000,1,...
                                        Info_enzyme,...
                                        Info_mRNA,...
                                        Info_protein,...
                                        Info_ribosome,...
                                        Info_tRNA);% kcat_glc*10000 used to eliminate glucose transporter which has been removed
            command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
            system(command,'-echo');
            fileName_out = 'Simulation.lp.out';
            [~,sol_status,sol_full] = ReadSoplexResult(fileName_out,model);

            if strcmp(sol_status,'optimal')
                [~,~,f_enzyme_inact] = CheckInactiveEnzyme(model,sol_full);
                disp(['f_enzyme_inact = ' num2str(f_enzyme_inact)]);
                tmp_protein = f_enzyme_inact - f_transporter * 0.46;
                if tmp_protein > 0
                    factor_high = factor_mid;
                else
                    factor_low = factor_mid;
                end
            else
                disp('No optimal solution.');
                factor_low = factor_mid;
            end
        end
        fileName = WriteLPSatFactor(model,mu,f,osenseStr,rxnID,factor_low,...
                                    f_transporter,kcat_glc*10000,1,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,sol_status_low,sol_full_low] = ReadSoplexResult(fileName_out,model);

        fileName = WriteLPSatFactor(model,mu,f,osenseStr,rxnID,factor_high,...
                                    f_transporter,kcat_glc*10000,1,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,sol_status_high,sol_full_high] = ReadSoplexResult(fileName_out,model);
        
        if strcmp(sol_status_low,'optimal')
            [~,~,f_enzyme_inact_low] = CheckInactiveEnzyme(model,sol_full_low);
            tmp_protein_low = f_enzyme_inact_low - f_transporter * 0.46;
            if tmp_protein_low < 0.0001 && tmp_protein_low > -0.0001
                factor_k_new = factor_low;
                sol_full_new = sol_full_low;
            else
                if strcmp(sol_status_high,'optimal')
                    [~,~,f_enzyme_inact_high] = CheckInactiveEnzyme(model,sol_full_high);
                    tmp_protein_high = f_enzyme_inact_high - f_transporter * 0.46;
                    if tmp_protein_high < 0.0001 && tmp_protein_high > -0.0001
                        factor_k_new = factor_high;
                        sol_full_new = sol_full_high;
                    else
                        factor_k_new = nan;
                        sol_full_new = zeros(length(model.rxns),1);
                    end
                else
                    factor_k_new = nan;
                    sol_full_new = zeros(length(model.rxns),1);
                end
            end
        else
            if strcmp(sol_status_high,'optimal')
                [~,~,f_enzyme_inact_high] = CheckInactiveEnzyme(model,sol_full_high);
                tmp_protein_high = f_enzyme_inact_high - f_transporter * 0.46;
                if tmp_protein_high < 0.0001 && tmp_protein_high > -0.0001
                    factor_k_new = factor_high;
                    sol_full_new = sol_full_high;
                else
                    factor_k_new = nan;
                    sol_full_new = zeros(length(model.rxns),1);
                end
            else
                factor_k_new = nan;
                sol_full_new = zeros(length(model.rxns),1);
            end
        end
    end

elseif Check == 1
    
    factor_org = 1;

    fileName = WriteLPSatFactor(model,mu,f,osenseStr,rxnID,factor_org,...
                                f_transporter,kcat_glc,factor_org,...
                                Info_enzyme,...
                                Info_mRNA,...
                                Info_protein,...
                                Info_ribosome,...
                                Info_tRNA);
    command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
    system(command,'-echo');
    fileName_out = 'Simulation.lp.out';
    [~,sol_status,sol_full] = ReadSoplexResult(fileName_out,model);

    if strcmp(sol_status,'optimal')
        [~,~,f_enzyme_inact] = CheckInactiveEnzyme(model,sol_full);
    else 
        error('No solution is found for the original simulation.');
    end

    tmp_protein = f_enzyme_inact;

    if tmp_protein == 0
        factor_k_new = factor_org;
        sol_full_new = sol_full;
    else
        factor_high = 1;
        factor_low = 0;
        while factor_high-factor_low > 0.0001 || (tmp_protein > 0 && factor_high-factor_low > 0.000001) || ~strcmp(sol_status,'optimal')
            factor_mid = (factor_low+factor_high)/2;
            disp(['kcat factor = ' num2str(factor_mid)]);

            fileName = WriteLPSatFactor(model,mu,f,osenseStr,rxnID,factor_mid,...
                                        f_transporter,kcat_glc,factor_mid,...
                                        Info_enzyme,...
                                        Info_mRNA,...
                                        Info_protein,...
                                        Info_ribosome,...
                                        Info_tRNA);
            command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
            system(command,'-echo');
            fileName_out = 'Simulation.lp.out';
            [~,sol_status,sol_full] = ReadSoplexResult(fileName_out,model);

            if strcmp(sol_status,'optimal')
                [~,~,f_enzyme_inact] = CheckInactiveEnzyme(model,sol_full);
                disp(['f_enzyme_inact = ' num2str(f_enzyme_inact)]);
                tmp_protein = f_enzyme_inact;
                if tmp_protein > 0
                    factor_high = factor_mid;
                else
                    factor_low = factor_mid;
                end
            else
                disp('No optimal solution.');
                factor_low = factor_mid;
            end
        end
        fileName = WriteLPSatFactor(model,mu,f,osenseStr,rxnID,factor_low,...
                                    f_transporter,kcat_glc,factor_low,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,sol_status_low,sol_full_low] = ReadSoplexResult(fileName_out,model);

        fileName = WriteLPSatFactor(model,mu,f,osenseStr,rxnID,factor_high,...
                                    f_transporter,kcat_glc,factor_high,...
                                    Info_enzyme,...
                                    Info_mRNA,...
                                    Info_protein,...
                                    Info_ribosome,...
                                    Info_tRNA);
        command = sprintf('/Users/cheyu/build/bin/soplex -s0 -g5 -f1e-18 -o1e-18 -x -q -c --int:readmode=1 --int:solvemode=2 --int:checkmode=2 %s > %s.out %s',fileName,fileName);
        system(command,'-echo');
        fileName_out = 'Simulation.lp.out';
        [~,sol_status_high,sol_full_high] = ReadSoplexResult(fileName_out,model);

        if strcmp(sol_status_low,'optimal')
            [~,~,f_enzyme_inact_low] = CheckInactiveEnzyme(model,sol_full_low);
            tmp_protein_low = f_enzyme_inact_low;
            if tmp_protein_low < 0.0001 && tmp_protein_low > -0.0001
                factor_k_new = factor_low;
                sol_full_new = sol_full_low;
            else
                if strcmp(sol_status_high,'optimal')
                    [~,~,f_enzyme_inact_high] = CheckInactiveEnzyme(model,sol_full_high);
                    tmp_protein_high = f_enzyme_inact_high - f_transporter * 0.46;
                    if tmp_protein_high < 0.0001 && tmp_protein_high > -0.0001
                        factor_k_new = factor_high;
                        sol_full_new = sol_full_high;
                    else
                        factor_k_new = nan;
                        sol_full_new = zeros(length(model.rxns),1);
                    end
                else
                    factor_k_new = nan;
                    sol_full_new = zeros(length(model.rxns),1);
                end
            end
        else
            if strcmp(sol_status_high,'optimal')
                [~,~,f_enzyme_inact_high] = CheckInactiveEnzyme(model,sol_full_high);
                tmp_protein_high = f_enzyme_inact_high;
                if tmp_protein_high < 0.0001 && tmp_protein_high > -0.0001
                    factor_k_new = factor_high;
                    sol_full_new = sol_full_high;
                else
                    factor_k_new = nan;
                    sol_full_new = zeros(length(model.rxns),1);
                end
            else
                factor_k_new = nan;
                sol_full_new = zeros(length(model.rxns),1);
            end
        end
    end
end

