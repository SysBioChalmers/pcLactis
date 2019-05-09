%% AddMetabolicRxn
%   This function will generate metabolic reactions.
function [M_Metabolism,...
          E_Enzyme_formation,...
          E_Enzyme_dilution,...
          E_Protein_degradation,...
          RxnGPR] = AddMetabolicRxn(Metabolic_Matrix,...
                                    Functional_protein,...
                                    E_Enzyme_formation,...
                                    E_Enzyme_dilution,...
                                    E_Protein_degradation)

% Find reactions with GPRs.
RxnGPR = struct();
RxnGPR.Rxn = {};
RxnGPR.GPR = {};
RxnGPR.newGPR = {};

Rxn_list = unique(Metabolic_Matrix.RxnList);
for i = 1:length(Rxn_list)
    Rxn_name = Rxn_list{i};
    [~, index] = ismember(Rxn_name,Metabolic_Matrix.RxnList);
    GPR = Metabolic_Matrix.CatalystList{index};
    if ~isempty(GPR)%if has GPR
        newGPR = {strcat(Rxn_name(3:end),'_Enzyme_c')};
        RxnGPR.Rxn = cat(1,RxnGPR.Rxn,Rxn_name);
        RxnGPR.GPR = cat(1,RxnGPR.GPR,GPR);
        RxnGPR.newGPR = cat(1,RxnGPR.newGPR,newGPR);
    end
end

% Update Metabolic_Matrix to M_Metabolism.
M_Metabolism = Metabolic_Matrix;
for i = 1:length(M_Metabolism.CatalystList)
    Cata = M_Metabolism.CatalystList{i};
    if ~isempty(Cata)%if has GPR
        Rxn_name = M_Metabolism.RxnList{i};
        M_Metabolism.CatalystList{i} = strcat(Rxn_name(3:end),'_Enzyme_c');
    end
end

% Add enzyme formation, dilution and protein degradation reactions.
for i = 1:length(RxnGPR.Rxn)
    Rxn_name = RxnGPR.Rxn{i};
    subGPR = RxnGPR.GPR{i};
    cplxGPR = RxnGPR.newGPR{i};

    if ~contains(subGPR,'and')%if has a single enzyme
            
        %add metabolic enzyme formation reactions    
        RxnID = strcat(Rxn_name,'_Enzyme');
        Rev = 0;
        Subproc = 'Enzyme formation';
        Catalyst = '';
        
        Sub = FindFunctionalProtein(subGPR,Functional_protein);
        Compo = {Sub;cplxGPR};
        Coeff = [-1;1];
        E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
        
        %add protein degradation reactions    
        RxnID = strcat(Rxn_name,'_Enzyme_degradation');
        Subproc = 'Protein degradation';
        Sub = strcat(Sub(1:end-2),'_degradation_c');
        Compo = {cplxGPR;Sub};
        Coeff = [-1;1];
        E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
                                 
    else%if has a complex
        %add metabolic enzyme formation reactions    
        RxnID = strcat(Rxn_name,'_Enzyme');
        Rev = 0;
        Subproc = 'Enzyme formation';
        Catalyst = '';
        
        Sub_raw = transpose(strsplit(subGPR,' and '));
        Sub = cellfun(@(x) FindFunctionalProtein(x,Functional_protein),...
                      Sub_raw,'UniformOutput',false);
        Compo = cat(1,Sub,{cplxGPR});
        n_sub = length(Sub);
        Coeff = [-1*ones(n_sub,1);1];
        E_Enzyme_formation = MatrixOneReactionAdding(E_Enzyme_formation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
                                 
        %add protein degradation reactions    
        RxnID = strcat(Rxn_name,'_Enzyme_degradation');
        Subproc = 'Protein degradation';
        Sub_new = cellfun(@(x) strcat(x(1:end-2),'_degradation_c'),...
                      Sub,'UniformOutput',false);
        Compo = cat(1,{cplxGPR},Sub_new);
        n_sub = length(Sub_new);
        Coeff = [-1;1*ones(n_sub,1)];
        E_Protein_degradation = MatrixOneReactionAdding(E_Protein_degradation,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
    end
    
    %add dilution reactions
    RxnID = strcat('R_dilution_',cplxGPR(1:end-2));
    Compo = {cplxGPR};
	Coeff = -1;
	Rev = 0;
	Subproc = 'Enzyme dilution';
	Catalyst = '';
	E_Enzyme_dilution = MatrixOneReactionAdding(E_Enzyme_dilution,...
                                     RxnID,Compo,Coeff,Rev,Subproc,Catalyst);
end

