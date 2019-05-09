function fileName = WriteLPforM(model)

fileName = sprintf('Simulation_M.lp');
fptr = fopen(fileName,'w');

% Set objective function (here is to minimize glucose exchange reaction)
fprintf(fptr,'Maximize\n');
index_obj = find(ismember(model.rxns,'R_M_EX_glc_LPAREN_e_RPAREN_'));
fprintf(fptr,'obj: X%d\n',index_obj);
fprintf(fptr,'Subject To\n');

for i = 1:numel(model.mets)
    j = find(full(model.S(i,:)));
    for m = 1:numel(j)
        s = full(model.S(i,j(m)));
        if mod(m,200) == 0
            sep = newline;
        else
            sep = '';
        end
        if m == 1
           eq = sprintf('%.15f X%d',s,j(m));
        else
           if s>0
               eq = sprintf('%s + %.15f X%d%c',eq,s,j(m),sep);
           else
               eq = sprintf('%s %.15f X%d%c',eq,s,j(m),sep);
           end
        end
    end
    fprintf(fptr,'C%d: %s = 0\n',i,eq);
end

fprintf(fptr,'Bounds\n');

for i = 1:length(model.rxns)
	if model.ub(i) >= 100
        fprintf(fptr,'%f <= X%d <= +infinity\n',model.lb(i),i);
    else
        fprintf(fptr,'%f <= X%d <= %f\n',model.lb(i),i,model.ub(i));
	end
end

fprintf(fptr,'End\n');
fclose(fptr);

command =sprintf('/Users/cheyu/build/bin/soplex  -s0 -x -q -c %s > %s.out %s',fileName,fileName);
system(command,'-echo');