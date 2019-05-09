%% This function is from GECKO.
 function org_index = find_inKEGG(org_name,names)
     org_index      = find(strcmpi(org_name,names));
     if isempty(org_index)
         i=1;
         while isempty(org_index) && i<length(names)
             str = names{i};
             if strcmpi(org_name(1:strfind(org_name,' ')-1),...
                 str(1:strfind(str,' ')-1))
                 org_index = i;
             end
             i = i+1;
         end
         if isempty(org_index);org_index = '*';end
     end
 end