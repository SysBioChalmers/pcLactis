%% AddrRNAModificationRxn
function E_rRNA_modification = AddrRNAModificationRxn

% Generate matrix.
E_rRNA_modification = MatrixGeneration;

% Import rRNA modification data.
cd Data;
[~, ~, rRNAmod] = xlsread('RNA_modification.xlsx','rRNA');
cd ../;

% Add rRNA modification reactions.
data = struct();
data.type = rRNAmod(2:end,1);
data.cata = rRNAmod(2:end,4);
data.rxn = rRNAmod(2:end,5);

unique_type = unique(data.type);

for i = 1:length(unique_type)
    tpye = unique_type{i};
    n_type = find(strcmp(data.type,tpye));
    n = length(n_type);
    for j = 1
        RxnID = strcat('R_rRNA_modification_rRNA_',tpye,'_modified_',num2str(j));
        rRNAsub = strcat('rRNA_',tpye,'_unmodified_c');
        rRNAprod = strcat('rRNA_',tpye,'_modified_',num2str(j),'_c');
        
        E_rRNA_modification = AddrRNAMdfRxn(E_rRNA_modification,...
                                            RxnID,rRNAsub,rRNAprod,...
                                            data,n_type,j);
    end
        
    for j = 2:n-1
        RxnID = strcat('R_rRNA_modification_rRNA_',tpye,'_modified_',num2str(j));
        rRNAsub = strcat('rRNA_',tpye,'_modified_',num2str(j-1),'_c');
        rRNAprod = strcat('rRNA_',tpye,'_modified_',num2str(j),'_c');
        
        E_rRNA_modification = AddrRNAMdfRxn(E_rRNA_modification,...
                                            RxnID,rRNAsub,rRNAprod,...
                                            data,n_type,j);
    end
        
	for j = n
        RxnID = strcat('R_rRNA_modification_rRNA_',tpye);
        rRNAsub = strcat('rRNA_',tpye,'_modified_',num2str(j-1),'_c');
        rRNAprod = strcat('rRNA_',tpye,'_c');
        
        E_rRNA_modification = AddrRNAMdfRxn(E_rRNA_modification,...
                                            RxnID,rRNAsub,rRNAprod,...
                                            data,n_type,j);

	end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_rRNA_modification = AddrRNAMdfRxn(E_rRNA_modification,...
                                             RxnID,rRNAsub,rRNAprod,...
                                             data,n_type,j)
if all(isnan(data.rxn{n_type(j)}))%no other metabolites involved
	if all(isnan(data.cata{n_type(j)}))%no catalyst involved
        Catalyst = '';
        E_rRNA_modification = AddRxnForNoMet(E_rRNA_modification,...
                                             RxnID,rRNAsub,rRNAprod,...
                                             Catalyst);
	else%catalyst involved
        %add catalyst formation reaction
        Cata = data.cata{n_type(j)};
        if contains(Cata,'and')
            Cata = strrep(Cata,' and ','_');
        end
        %add rRNA modification reaction
        Catalyst = strcat('RNAmod_',Cata,'_Enzyme_c');
        E_rRNA_modification = AddRxnForNoMet(E_rRNA_modification,...
                                             RxnID,rRNAsub,rRNAprod,...
                                             Catalyst);
	end
else%other metabolites involved
	MetInfo = extractBetween(data.rxn{n_type(j)},'(',')');
	Metname = MetInfo{1};
    Metname = strrep(Metname,', ',' ');
    Metname = transpose(strsplit(Metname));    
    Metcoeff = transpose(str2num(MetInfo{2}));
            
	if all(isnan(data.cata{n_type(j)}))%no catalyst involved
        Catalyst = '';
        E_rRNA_modification = AddRxnForMet(E_rRNA_modification,...
                                           RxnID,rRNAsub,rRNAprod,...
                                           Metname,Metcoeff,Catalyst);
	else%catalyst involved
        %add catalyst formation reaction
        Cata = data.cata{n_type(j)};
        if contains(Cata,'and')
            Cata = strrep(Cata,' and ','_');
        end
        %add tRNA modification reaction
        Catalyst = strcat('RNAmod_',Cata,'_Enzyme_c');
        E_rRNA_modification = AddRxnForMet(E_rRNA_modification,...
                                           RxnID,rRNAsub,rRNAprod,...
                                           Metname,Metcoeff,Catalyst);
	end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_rRNA_modification = AddRxnForNoMet(E_rRNA_modification,...
                                              RxnID,rRNAsub,rRNAprod,...
                                              Catalyst)
Compo = {rRNAsub;rRNAprod};
Coeff = [-1;1];
Rev = 0;
Subproc = 'rRNA modification';

E_rRNA_modification = MatrixOneReactionAdding(E_rRNA_modification,...
                                      RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function E_rRNA_modification = AddRxnForMet(E_rRNA_modification,...
                                            RxnID,rRNAsub,rRNAprod,...
                                            Metname,Metcoeff,Catalyst)
Compo = cat(1,{rRNAsub;rRNAprod},Metname);
Coeff = [-1;1;Metcoeff];
Rev = 0;
Subproc = 'rRNA modification';

E_rRNA_modification = MatrixOneReactionAdding(E_rRNA_modification,...
                                      RxnID,Compo,Coeff,Rev,Subproc,Catalyst);

end
