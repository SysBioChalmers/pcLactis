function writeClusterFile

load('param_list.mat');
kmax = length(param_list);
k = 1:420:kmax;
for i = 1:length(k)
    subfileName = ['sub_up_1000_',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH --time 0-12:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=cheyu@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB intel/2018b GMP\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    for m = 1:20
        if i+21*(m-1) < kmax
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(21*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(21*m-1) '\n']);
        end
    end
    for m = 1:20
        if i+21*(m-1) < kmax
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "kcat_sensitivity(2,1000,$a',num2str(m),',$b',num2str(m),')" &\n']);
            fprintf(fptr,'sleep 60s\n');
        end
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end

for i = 1:length(k)
    subfileName = ['sub_up_10_',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH --time 0-12:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=cheyu@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB intel/2018b GMP\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    for m = 1:20
        if i+21*(m-1) < kmax
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(21*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(21*m-1) '\n']);
        end
    end
    for m = 1:20
        if i+21*(m-1) < kmax
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "kcat_sensitivity(2,10,$a',num2str(m),',$b',num2str(m),')" &\n']);
            fprintf(fptr,'sleep 60s\n');
        end
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end

for i = 1:length(k)
    subfileName = ['sub_dn_1000_',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH --time 0-12:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=cheyu@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB intel/2018b GMP\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    for m = 1:20
        if i+21*(m-1) < kmax
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(21*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(21*m-1) '\n']);
        end
    end
    for m = 1:20
        if i+21*(m-1) < kmax
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "kcat_sensitivity(0.5,1000,$a',num2str(m),',$b',num2str(m),')" &\n']);
            fprintf(fptr,'sleep 60s\n');
        end
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end

for i = 1:length(k)
    subfileName = ['sub_dn_10_',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH --time 0-12:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=cheyu@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB intel/2018b GMP\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    for m = 1:20
        if i+21*(m-1) < kmax
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(21*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(21*m-1) '\n']);
        end
    end
    for m = 1:20
        if i+21*(m-1) < kmax
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "kcat_sensitivity(0.5,10,$a',num2str(m),',$b',num2str(m),')" &\n']);
            fprintf(fptr,'sleep 60s\n');
        end
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end

for i = 1:length(k)
    subfileName = ['sub_dn_1747_',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH --time 0-12:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=cheyu@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB intel/2018b GMP\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    for m = 1:20
        if i+21*(m-1) < kmax
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(21*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(21*m-1) '\n']);
        end
    end
    for m = 1:20
        if i+21*(m-1) < kmax
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "kcat_sensitivity(0.5,17.47,$a',num2str(m),',$b',num2str(m),')" &\n']);
            fprintf(fptr,'sleep 60s\n');
        end
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end

for i = 1:length(k)
    subfileName = ['sub_up_1747_',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH --time 0-12:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=cheyu@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB intel/2018b GMP\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    for m = 1:20
        if i+21*(m-1) < kmax
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(21*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(21*m-1) '\n']);
        end
    end
    for m = 1:20
        if i+21*(m-1) < kmax
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "kcat_sensitivity(2,17.47,$a',num2str(m),',$b',num2str(m),')" &\n']);
            fprintf(fptr,'sleep 60s\n');
        end
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end