%% This function was from GECKO but modified a little. 
 
 function [dir,tot] =iterativeMatch(EC,MW,subs,i,j,KCATcell,dir,tot,m_matrix,...
                                    name,phylDist,org_index,SAcell,Rxns)
 %Will iteratively try to match the EC number to some registry in BRENDA,
 %using each time one additional wildcard.
 
 EC      = strsplit(EC,' ');
 kcat    = zeros(size(EC));
 origin  = zeros(size(EC));
 matches = zeros(size(EC));
 wc_num  = ones(size(EC)).*1000;
 for k = 1:length(EC)
     success  = false;
     while ~success
         %Atempt match:
         [kcat(k),origin(k),matches(k)] = mainMatch(EC{k},MW,subs,KCATcell,...
                                                    m_matrix,i,name,phylDist,...
                                                    org_index,SAcell,Rxns);
         %If any match found, ends. If not, introduces one extra wild card and
         %tries again:
         if origin(k) > 0
             success   = true;
             wc_num(k) = sum(EC{k}=='-');
         else
             dot_pos  = [2 strfind(EC{k},'.')];
             wild_num = sum(EC{k}=='-');
             wc_text  = '-.-.-.-';
             EC{k}    = [EC{k}(1:dot_pos(4-wild_num)) wc_text(1:2*wild_num+1)];
         end
     end
 end
 
 if sum(origin) > 0
     %For more than one EC: Choose the maximum value among the ones with the
     %less amount of wildcards and the better origin:
     best_pos   = (wc_num == min(wc_num));
     new_origin = origin(best_pos);
     best_pos   = (origin == min(new_origin(new_origin~=0)));
     max_pos    = find(kcat == max(kcat(best_pos)));
     wc_num     = wc_num(max_pos(1));
     origin     = origin(max_pos(1));
     matches    = matches(max_pos(1));
     kcat       = kcat(max_pos(1));
     
     %Update dir and tot:
     dir.kcats(i,j)   = kcat;
     dir.org_s(i,j)   = matches*(origin == 1);
     dir.rest_s(i,j)  = matches*(origin == 2);
     dir.org_ns(i,j)  = matches*(origin == 3);
     dir.org_sa(i,j)  = matches*(origin == 4);
     dir.rest_ns(i,j) = matches*(origin == 5);    
     dir.rest_sa(i,j) = matches*(origin == 6);
     tot.org_s        = tot.org_s   + (origin == 1);
     tot.rest_s       = tot.rest_s  + (origin == 2);
     tot.org_ns       = tot.org_ns  + (origin == 3);
     tot.org_sa       = tot.org_sa  + (origin == 4);
     tot.rest_ns      = tot.rest_ns + (origin == 5);    
     tot.rest_sa      = tot.rest_sa + (origin == 6);
     tot.wc0          = tot.wc0     + (wc_num == 0);
     tot.wc1          = tot.wc1     + (wc_num == 1);
     tot.wc2          = tot.wc2     + (wc_num == 2);
     tot.wc3          = tot.wc3     + (wc_num == 3);
     tot.wc4          = tot.wc4     + (wc_num == 4);
     tot.queries      = tot.queries + 1;
     tot.matrix(origin,wc_num+1) = tot.matrix(origin,wc_num+1) + 1;
 end
 
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 function [kcat,origin,matches] = mainMatch(EC,MW,subs,KCATcell,m_matrix,i,...
                                      name,phylDist,org_index,SAcell,Rxns)
                                                                   
 % Matching function prioritizing organism and substrate specificity when 
 % available.
 
 origin = 0;
 %First try to match organism and substrate:
 [kcat,matches] = matchKcat(EC,MW,subs,KCATcell,name,true,false,m_matrix,i,...
                            phylDist,org_index,SAcell,Rxns);                      
 if matches > 0
     origin = 1;
 %If no match, try the closest organism but match the substrate:
 else   
    [kcat,matches] = matchKcat(EC,MW,subs,KCATcell,'',true,false,m_matrix,i,...
                               phylDist,org_index,SAcell,Rxns);
     if matches > 0
         origin = 2;
     %If no match, try to match organism but with any substrate:
     else
         [kcat,matches] = matchKcat(EC,MW,subs,KCATcell,name,false,false,...
                                    m_matrix,i,phylDist,org_index,SAcell,Rxns);
         if matches > 0
             origin = 3;
         %If no match, try to match organism but for any substrate (SA*MW):
         else
              [kcat,matches] = matchKcat(EC,MW,subs,KCATcell,name,false,...
                                         true,m_matrix,i,phylDist,org_index,...
                                         SAcell,Rxns);
              if matches > 0
                  origin = 4; 
             %If no match, try any organism and any substrate:
              else
                 [kcat,matches] = matchKcat(EC,MW,subs,KCATcell,'',false,...
                                            false,m_matrix,i,phylDist,...
                                            org_index,SAcell,Rxns);
                 if matches > 0
                     origin = 5;
                 %Again if no match, look for any org and SA*MW    
                  else
                      [kcat,matches] = matchKcat(EC,MW,subs,KCATcell,'',...
                                                 false,true,m_matrix,i,phylDist,...
                                                 org_index, SAcell,Rxns);
                      if matches > 0
                          origin = 6;
                      end
                 end
                         
              end    
         end
     end
 end
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [kcat,matches] = matchKcat(EC,MW,subs,KCATcell,organism,...
                                     substrate,SA,m_matrix,i,phylDist,...
                                                         org_index,SAcell,Rxns)
                  
 %Will go through BRENDA and will record any match. Afterwards, it will
 %return the average value and the number of matches attained.
 kcat    = [];
 matches = 0;
 %Relaxes matching if wild cards are present:
 wild     = false;
 wild_pos = strfind(EC,'-');
 if ~isempty(wild_pos)
     EC   = EC(1:wild_pos(1)-1);
     wild = true;
 end

 if SA
     EC_indexes = extract_indexes(EC,SAcell{1},[],SAcell{2},subs,substrate,...
                                  organism,org_index,phylDist,wild);
     kcat       = SAcell{3}(EC_indexes);
     org_cell   = SAcell{2}(EC_indexes);
     MW_BRENDA  = SAcell{4}(EC_indexes);
 else
     EC_indexes = extract_indexes(EC,KCATcell{1},KCATcell{2},KCATcell{3},...
                                  subs,substrate,organism,org_index,...
                                                      phylDist,wild);
     if substrate
         for j = 1:length(EC_indexes)
             indx = EC_indexes(j);
             for k = 1:length(subs)
                 l = boolean(strcmpi(m_matrix.MetName,subs{k}).*(strcmp(Rxns{i},m_matrix.RxnList)));
                 if ~isempty(subs{k}) && strcmpi(subs{k},KCATcell{2}(indx))
                     if KCATcell{4}(indx) > 0 
                         coeff = min(abs(m_matrix.CoeffList(l,1)));
                         kcat  = [kcat;KCATcell{4}(indx)/coeff];
                     end
                 end
             end
         end
     else
         kcat = KCATcell{4}(EC_indexes);
     end
 end                         
 %Return maximum value:
 if isempty(kcat)
     kcat = 0;
 else
     matches        = length(kcat);
     [kcat,MaxIndx] = max(kcat);
%      if SA
%          % If the match correspond to a SA*Mw value for the model's
%          % organism the kcat will be corrected with the sequence based Mw
%          if strcmpi(organism,org_cell(MaxIndx))
%              kcat = kcat*MW/MW_BRENDA(MaxIndx);
%          end
%      end        
 end
 %Avoid SA*Mw values over the diffusion limit rate  [Bar-Even et al. 2011]
 if kcat>(1E7*3600)
     kcat = 1E7*3600;
 end
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Extract the indexes of the entries in the BRENDA data that meet the 
 %conditions specified by the search criteria
 function EC_indexes = extract_indexes(EC,EC_cell,subs_cell,orgs_cell,subs,...
                                       substrate,organism, org_index,...
                                       phylDist,wild)
 
 EC_indexes = [];
 if wild
   for j=1:length(EC_cell)
        if strfind(EC_cell{j},EC)==1
            EC_indexes = [EC_indexes,j];
        end
    end   
 else
    EC_indexes = transpose(find(strcmpi(EC,EC_cell)));  
 end
 %If substrate=true then it will extract only the substrates appereances 
 %indexes in the EC subset from the BRENDA cell array
 if substrate
     Subs_indexes = [];
     for l = 1:length(subs)
         if ~isempty(subs(l))
             Subs_indexes = horzcat(Subs_indexes,EC_indexes(strcmpi(subs(l),...
                                    subs_cell(EC_indexes))));          
         end
     end
     EC_indexes = Subs_indexes;    
 end
 
 EC_orgs = orgs_cell(EC_indexes);
 %If specific organism values are requested looks for all the organism
 %repetitions on the subset BRENDA cell array(EC_indexes)
 if string(organism) ~= ''  
     EC_indexes = EC_indexes(strcmpi(string(organism),EC_orgs));
 
 %If KEGG code was assigned to the organism (model) then it will look for   
 %the Kcat value for the closest organism
 elseif org_index~='*' %&& org_index~=''
     KEGG_indexes = [];temp = [];
     
     %For relating a phyl dist between the modelled organism and the organisms
     %on the BRENDA cell array it should search for a KEGG code for each of 
     %these 
     for j=1:length(EC_indexes)
 
         %Assigns a KEGG index for those found on the KEGG struct
         orgs_index = find(strcmpi(orgs_cell(EC_indexes(j)),phylDist.names),1);
         if ~isempty(orgs_index)
             KEGG_indexes = [KEGG_indexes; orgs_index];
             temp         = [temp;EC_indexes(j)];
         %For values related to organisms without KEGG code, then it will
         %look for KEGG code for the first organism with the same genus
         else
             k=1;
             while isempty(orgs_index) && k<length(phylDist.names)
                 str = phylDist.names{k};
                 org = orgs_cell{EC_indexes(j)};
                 if strcmpi(org(1:strfind(org,' ')-1),str(1:strfind(str,' ')-1))
                     orgs_index   = k; 
                     KEGG_indexes = [KEGG_indexes;k];
                     temp         = [temp;EC_indexes(j)];
                 end
                 k = k+1;
             end
         end
     end
     %Update the EC_indexes cell array
     EC_indexes = temp;
     %Looks for the taxonomically closest organism and saves the index of
     %its appearences in the BRENDA cell
     if ~isempty(EC_indexes)
         distances  = num2cell(phylDist.distMat(org_index,:));
         distances  = distances(KEGG_indexes);
         EC_indexes = EC_indexes(cell2mat(distances) ==...
                                                 min(cell2mat(distances)));                     
     end
 end
 
 end 


