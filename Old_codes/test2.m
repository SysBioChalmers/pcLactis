
yatp_1 = 3;
yatp_2 = 2;
% yatp_3 = 3;
yatp_3 = 1.67;
protein_efficiency_1 = 456.073;
protein_efficiency_2 = 896.77;
protein_efficiency_3 = 230; % 230
totprotcost_1 = yatp_1/protein_efficiency_1;
totprotcost_2 = yatp_2/protein_efficiency_2;
totprotcost_3 = yatp_3/protein_efficiency_3;


%% Small model simulation

% Estimate ATP consumption for growth, which is the sum of ATP for GAM,
% R_M_PROTS_LLA_v3, R_M_RNAS_LLA, and R_M_DNAS_LLA. The others can be
% neglected. 
GAM_MModel = 33.783;
newGAM = GAM_MModel + 18.0895 + 0.1316 + 0.10138;
NGAM_MModel = 3.72;

tot_proteome = 0.46*0.1;


% Estimate upper limit on arg uptake. The small model does not account for
% biomass formation, so arginine is only used to produce ATP.
argUB = 1.6568 - 0.1764; % total arg consumption minus arg composition in biomass.

% Simulations
mu_list = 0.01:0.01:0.63;
flux = zeros(2,length(mu_list));
unused_prot = zeros(1,length(mu_list));
for i = 1:length(mu_list)
    mu = mu_list(i);

    A = [-yatp_1 -yatp_2;totprotcost_1 totprotcost_2];
    
    b = [-newGAM*mu-NGAM_MModel+yatp_3*argUB*mu tot_proteome-totprotcost_3*argUB*mu];

    Aeq = [];
    beq = [];

    lb = [0 0];
    ub = [100 100];

    f = [1 1];

    x = linprog(f,A,b,Aeq,beq,lb,ub);
    
    flux(:,i) = x;
    used_prot = totprotcost_1*x(1)+totprotcost_2*x(2)+totprotcost_3*argUB*mu;
    unused_prot(i) = tot_proteome - round(used_prot,6);
end

clear i mu A b Aeq beq lb ub f x factor used_prot;

%% Sensitivity analysis

% increasefactor = 0.01;
% 
% % 1. reference
% res_ref = zeros(4,length(flux)); %[X;Y;Z;ref_mu]
% 
% A = [totprotcost_1 totprotcost_2 totprotcost_3];
% b = tot_proteome;
% Aeq = [];
% beq = [];
% f = [-yatp_1 -yatp_2 -yatp_3];
% for i = 1:length(flux)
%     lb = flux(:,i)';
%     ub = flux(:,i)';
%     [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub);
%     res_ref(1:3,i) = x;
%     res_ref(4,i) = -(fval+NGAM_MModel)/newGAM;
% end
% 
% 
% % 2. increase bound of arg uptake
% res_arg = zeros(5,length(res_ref)); %[ref_arg;new_arg;ref_mu;new_mu;qs]
% 
% f = [-yatp_1 -yatp_2 -yatp_3];
% for i = 1:length(res_ref)
%     A = [1 1 0;totprotcost_1 totprotcost_2 totprotcost_3];
%     b = [sum(res_ref(1:2,i)) tot_proteome];
%     Aeq = [];
%     beq = [];
% 	lb = [0 0 0];
%     ub = [100 100 argUB*res_ref(4,i)*(1+increasefactor)];
%     [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub);
%     res_arg(1,i) = argUB*res_ref(4,i);
%     res_arg(2,i) = argUB*res_ref(4,i)*(1+increasefactor);
%     res_arg(4,i) = -(fval+NGAM_MModel)/newGAM;
%     res_arg(5,i) = x(1)+x(2);
% end
% res_arg(3,:) = res_ref(4,:);
% sensitivity_arg = (res_arg(4,:)-res_arg(3,:))./(res_arg(2,:)-res_arg(1,:));
% sensitivity_arg = round(sensitivity_arg,6);
% sensitivity_arg = sensitivity_arg.*round((res_arg(1,:)./res_arg(3,:)),6);
% 
% % 3. increase bound of glucose uptake, similar to increase glucose
% % transporter.
% res_glc = zeros(5,length(res_ref)); %[ref_glc;new_glc;ref_mu;new_mu;qs]
% 
% f = [-yatp_1 -yatp_2 -yatp_3];
% for i = 1:length(res_ref)
%     A = [1 1 0;totprotcost_1 totprotcost_2 totprotcost_3];
%     b = [sum(res_ref(1:2,i))*(1+increasefactor) tot_proteome];
%     Aeq = [];
%     beq = [];
% 	lb = [0 0 0];
%     ub = [100 100 argUB*res_ref(4,i)];
%     [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub);
%     res_glc(1,i) = sum(res_ref(1:2,i));
%     res_glc(2,i) = sum(res_ref(1:2,i))*(1+increasefactor);
%     res_glc(4,i) = -(fval+NGAM_MModel)/newGAM;
%     res_glc(5,i) = sum(res_ref(1:2,i));
% end
% res_glc(3,:) = res_ref(4,:);
% sensitivity_glc = (res_glc(4,:)-res_glc(3,:))./(res_glc(2,:)-res_glc(1,:));
% sensitivity_glc = round(sensitivity_glc,6);
% sensitivity_glc = sensitivity_glc.*round((res_glc(1,:)./res_glc(3,:)),6);
% 
% % 4. increase constraint on proteome
% res_proteome = zeros(5,length(res_ref)); %[ref_proteome;new_proteome;ref_mu;new_mu;qs]
% 
% f = [-yatp_1 -yatp_2 -yatp_3];
% for i = 1:length(res_ref)
%     A = [1 1 0;totprotcost_1 totprotcost_2 totprotcost_3];
%     b = [sum(res_ref(1:2,i)) tot_proteome*(1+increasefactor)];
%     Aeq = [];
%     beq = [];
% 	lb = [0 0 0];
%     ub = [100 100 argUB*res_ref(4,i)];
%     [x, fval] = linprog(f,A,b,Aeq,beq,lb,ub);
%     res_proteome(1,i) = tot_proteome;
%     res_proteome(2,i) = tot_proteome*(1+increasefactor);
%     res_proteome(4,i) = -(fval+NGAM_MModel)/newGAM;
%     res_proteome(5,i) = x(1)+x(2);
% end
% res_proteome(3,:) = res_ref(4,:);
% sensitivity_proteome = (res_proteome(4,:)-res_proteome(3,:))./(res_proteome(2,:)-res_proteome(1,:));
% sensitivity_proteome = round(sensitivity_proteome,6);
% sensitivity_proteome = sensitivity_proteome.*round((res_proteome(1,:)./res_proteome(3,:)),6);





