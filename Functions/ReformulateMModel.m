%% ReformulateMModel 
function [Metabolic_Matrix, M_Model_total] = ReformulateMModel(M_Matrix_original)

% Generate a matrix for splitting isozymes.
M_Matrix_1 = MatrixGeneration;

% Assuming there are the following types of GPRs in the M model.
% 'geneA'; '(geneA or geneB or ...)'; '(geneA and geneB and ...)';
% '((geneA or geneB or ...) and (geneC or geneD or ...) and ...)';
% '((geneA and geneB and ...) or (geneC and geneD and ...) or ...)'.
% And no GPRs with more complicated structure,
% e.g., '((geneA or (geneB and geneC)) and (geneD or geneE))'.
% No such complicated GPRs in L. lactis M model, but if there are such
% GPRs, the code below should be modified.

for i = 1:length(M_Matrix_original.RxnList)
    z = strtrim(M_Matrix_original.CatalystList{i});%get GPR and delete leading and trailing whitespace
    x = length(find(cellfun('isempty',M_Matrix_1.RxnList) == 0));%count rows in M_Matrix_1
    if isempty(z)%if is spontaneous reaction
        M_Matrix_1.RxnList(x+1,1) = M_Matrix_original.RxnList(i,1);
        M_Matrix_1.CompoList(x+1,1) = M_Matrix_original.CompoList(i,1);
        M_Matrix_1.CoeffList(x+1,1) = M_Matrix_original.CoeffList(i,1);
        M_Matrix_1.CatalystList(x+1,1) = {z};
        M_Matrix_1.RevList(x+1,1) = M_Matrix_original.RevList(i,1);
        M_Matrix_1.SubprocList(x+1,1) = M_Matrix_original.SubprocList(i,1);
    else%if is enzymatic reaction
        if ~contains(z,'or') && ~contains(z,'and')%if has single enzyme
            z = strrep(z,'(',''); z = strrep(z,')',''); z = strtrim(z);
            %delete parentheses and leading and trailing whitespace if it has
            M_Matrix_1.RxnList(x+1,1) = M_Matrix_original.RxnList(i,1);
            M_Matrix_1.CompoList(x+1,1) = M_Matrix_original.CompoList(i,1);
            M_Matrix_1.CoeffList(x+1,1) = M_Matrix_original.CoeffList(i,1);
            M_Matrix_1.CatalystList(x+1,1) = {z};
            M_Matrix_1.RevList(x+1,1) = M_Matrix_original.RevList(i,1);
            M_Matrix_1.SubprocList(x+1,1) = M_Matrix_original.SubprocList(i,1);
        else%if has multiple genes, either isozymes or complexes
            if ~contains(z,'or')%if only has complex
                z = strrep(z,'(',''); z = strrep(z,')',''); z = strtrim(z);
                %delete parentheses and leading and trailing whitespace
                M_Matrix_1.RxnList(x+1,1) = M_Matrix_original.RxnList(i,1);
                M_Matrix_1.CompoList(x+1,1) = M_Matrix_original.CompoList(i,1);
                M_Matrix_1.CoeffList(x+1,1) = M_Matrix_original.CoeffList(i,1);
                M_Matrix_1.CatalystList(x+1,1) = {z};
                M_Matrix_1.RevList(x+1,1) = M_Matrix_original.RevList(i,1);
                M_Matrix_1.SubprocList(x+1,1) = M_Matrix_original.SubprocList(i,1);
            else%if has isozymes, and has or does not have complex
                if ~contains(z,'and')%if only has isozymes
                    z = strrep(z,'(',''); z = strrep(z,')','');
                    z = strrep(z,'or ',' '); z = strtrim(z);
                    z = strsplit(z);%split isozymes
                    for j = 1:length(z)%change reaction ID
                        M_Matrix_1.RxnList(x+j,1) = {strcat(M_Matrix_original.RxnList{i},'_',num2str(j))};
                    end
                    M_Matrix_1.CompoList(x+1:x+length(z),1) = M_Matrix_original.CompoList(i,1);
                    M_Matrix_1.CoeffList(x+1:x+length(z),1) = M_Matrix_original.CoeffList(i,1);
                    M_Matrix_1.CatalystList(x+1:x+length(z),1) = transpose(z);
                    M_Matrix_1.RevList(x+1:x+length(z),1) = M_Matrix_original.RevList(i,1);
                    M_Matrix_1.SubprocList(x+1:x+length(z),1) = M_Matrix_original.SubprocList(i,1);
                else%if has both isozymes and complex
                    if contains(z,') or')||contains(z,'or (')
                    % if GPR is made up of iso-complexes, e.g., '((cplxA) or (cplxB))'
                        z = strrep(z,') or',';'); z = strrep(z,'or (',';');
                        z = split(z,';');%split iso-complexes
                        z = strrep(z,'(',''); z = strrep(z,')',''); z = strtrim(z);
                        if contains(z,'or')
                            error('There are complicated GPRs in M model.')
                        else
                            for j = 1:length(z)%change reaction ID
                                M_Matrix_1.RxnList(x+j,1) = {strcat(M_Matrix_original.RxnList{i},'_',num2str(j))};
                            end
                            M_Matrix_1.CompoList(x+1:x+length(z),1) = M_Matrix_original.CompoList(i,1);
                            M_Matrix_1.CoeffList(x+1:x+length(z),1) = M_Matrix_original.CoeffList(i,1);
                            M_Matrix_1.CatalystList(x+1:x+length(z),1) = z;
                            M_Matrix_1.RevList(x+1:x+length(z),1) = M_Matrix_original.RevList(i,1);
                            M_Matrix_1.SubprocList(x+1:x+length(z),1) = M_Matrix_original.SubprocList(i,1);
                        end
                    elseif contains(z,') and')||contains(z,'and (')
                    % if GPR is a complex with isozymes, e.g., '((geneA or geneB) and (geneC or geneD))'
                        while z(1) == '('
                            z = z(2:end-1); z = strtrim(z);
                        end
                        z = {z};
                        while contains(z,' or ')
                        %The loop splits isozymes at a time 
                        %(e.g., (geneA or geneB) and next round for (geneC or geneD)) till no isozymes
                            isos = extractBetween(z,'(',')');%extract isozymes
                            iso = split(isos{1},' or ');%split isozymes
                            z_cell = cell(0,1);
                            for niso = 1:length(iso)
                                z_cell(niso,1) = {strrep(z,['(' isos{1} ')'],iso{niso})};
                            end
                            if ischar(z_cell{1})
                                z = z_cell;
                            else
                                z = z_cell{1};
                                for z_cell_raw = 2:length(z_cell)
                                    z = [z;z_cell{z_cell_raw}];
                                end
                            end
                        end
                        z = strrep(z,'  ',' '); z = strtrim(z);
                        for k = 1:length(z)%change reaction ID
                                M_Matrix_1.RxnList(x+k,1) = {strcat(M_Matrix_original.RxnList{i},'_',num2str(k))};
                        end
                        M_Matrix_1.CompoList(x+1:x+length(z),1) = M_Matrix_original.CompoList(i,1);
                        M_Matrix_1.CoeffList(x+1:x+length(z),1) = M_Matrix_original.CoeffList(i,1);
                        M_Matrix_1.CatalystList(x+1:x+length(z),1) = z;
                        M_Matrix_1.RevList(x+1:x+length(z),1) = M_Matrix_original.RevList(i,1);
                        M_Matrix_1.SubprocList(x+1:x+length(z),1) = M_Matrix_original.SubprocList(i,1);
                    else
                        error('There are complicated GPRs in M model.')
                    end
                end
            end
        end
    end
end

% Generate a matrix for splitting reversible reactions.
Metabolic_Matrix = MatrixGeneration;

for i = 1:length(M_Matrix_1.RevList)
    x = length(find(cellfun('isempty',Metabolic_Matrix.RxnList) == 0));%count rows in Metabolic_Matrix
    if isempty(M_Matrix_1.CatalystList{i})%if is spontaneous reaction
        Metabolic_Matrix.RxnList(x+1,1) = M_Matrix_1.RxnList(i,1);
        Metabolic_Matrix.CompoList(x+1,1) = M_Matrix_1.CompoList(i,1);
        Metabolic_Matrix.CoeffList(x+1,1) = M_Matrix_1.CoeffList(i,1);
        Metabolic_Matrix.CatalystList(x+1,1) = M_Matrix_1.CatalystList(i,1);
        Metabolic_Matrix.RevList(x+1,1) = M_Matrix_1.RevList(i,1);
        Metabolic_Matrix.SubprocList(x+1,1) = M_Matrix_1.SubprocList(i,1);
    else%if is enzymatic reaction
        if M_Matrix_1.RevList(i,1) == 0%if is not reversible
            Metabolic_Matrix.RxnList(x+1,1) = M_Matrix_1.RxnList(i,1);
            Metabolic_Matrix.CompoList(x+1,1) = M_Matrix_1.CompoList(i,1);
            Metabolic_Matrix.CoeffList(x+1,1) = M_Matrix_1.CoeffList(i,1);
            Metabolic_Matrix.CatalystList(x+1,1) = M_Matrix_1.CatalystList(i,1);
            Metabolic_Matrix.RevList(x+1,1) = M_Matrix_1.RevList(i,1);
            Metabolic_Matrix.SubprocList(x+1,1) = M_Matrix_1.SubprocList(i,1);
        else%if is reversible
            %add forward rows
            Metabolic_Matrix.RxnList(x+1,1) = {strcat(M_Matrix_1.RxnList{i},'_fwd')};%change reaction ID
            Metabolic_Matrix.CompoList(x+1,1) = M_Matrix_1.CompoList(i,1);
            Metabolic_Matrix.CoeffList(x+1,1) = M_Matrix_1.CoeffList(i,1);
            Metabolic_Matrix.CatalystList(x+1,1) = M_Matrix_1.CatalystList(i,1);
            Metabolic_Matrix.RevList(x+1,1) = 0;
            Metabolic_Matrix.SubprocList(x+1,1) = M_Matrix_1.SubprocList(i,1);
            %add reverse rows
            Metabolic_Matrix.RxnList(x+2,1) = {strcat(M_Matrix_1.RxnList{i},'_rvs')};%change reaction ID
            Metabolic_Matrix.CompoList(x+2,1) = M_Matrix_1.CompoList(i,1);
            Metabolic_Matrix.CoeffList(x+2,1) = -1*M_Matrix_1.CoeffList(i,1);
            Metabolic_Matrix.CatalystList(x+2,1) = M_Matrix_1.CatalystList(i,1);
            Metabolic_Matrix.RevList(x+2,1) = 0;
            Metabolic_Matrix.SubprocList(x+2,1) = M_Matrix_1.SubprocList(i,1);
        end
    end
end

M_Model_total = Matrix2CBModel(Metabolic_Matrix);


