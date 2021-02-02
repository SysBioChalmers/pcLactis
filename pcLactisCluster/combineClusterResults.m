function combineClusterResults
load('param_list.mat');
kmax = length(param_list);

cd ClusterResults/;

mu_list_kcat_sens_up_10 = zeros(0,0);
mu_list_kcat_sens_up_1000 = zeros(0,0);
mu_list_kcat_sens_dn_10 = zeros(0,0);
mu_list_kcat_sens_dn_1000 = zeros(0,0);

mu_list_kcat_sens_up_1747 = zeros(0,0);
mu_list_kcat_sens_dn_1747 = zeros(0,0);

k = 1:21:kmax;
for i = 1:length(k)
    display([num2str(i),'/',num2str(length(k))]);
    if i < length(k)
        file = ['kcat_sens_up_10_',num2str(k(i)),'_',num2str(k(i)+20),'.mat'];
        load(file);
        mu_list_kcat_sens_up_10 = [mu_list_kcat_sens_up_10;mu_list];

        file = ['kcat_sens_up_1000_',num2str(k(i)),'_',num2str(k(i)+20),'.mat'];
        load(file);
        mu_list_kcat_sens_up_1000 = [mu_list_kcat_sens_up_1000;mu_list];

        file = ['kcat_sens_dn_10_',num2str(k(i)),'_',num2str(k(i)+20),'.mat'];
        load(file);
        mu_list_kcat_sens_dn_10 = [mu_list_kcat_sens_dn_10;mu_list];

        file = ['kcat_sens_dn_1000_',num2str(k(i)),'_',num2str(k(i)+20),'.mat'];
        load(file);
        mu_list_kcat_sens_dn_1000 = [mu_list_kcat_sens_dn_1000;mu_list];
        
        file = ['kcat_sens_up_17.47_',num2str(k(i)),'_',num2str(k(i)+20),'.mat'];
        load(file);
        mu_list_kcat_sens_up_1747 = [mu_list_kcat_sens_up_1747;mu_list];
        
        file = ['kcat_sens_dn_17.47_',num2str(k(i)),'_',num2str(k(i)+20),'.mat'];
        load(file);
        mu_list_kcat_sens_dn_1747 = [mu_list_kcat_sens_dn_1747;mu_list];        
    elseif i == length(k)
        file = ['kcat_sens_up_10_',num2str(k(i)),'_',num2str(kmax),'.mat'];
        load(file);
        mu_list_kcat_sens_up_10 = [mu_list_kcat_sens_up_10;mu_list];

        file = ['kcat_sens_up_1000_',num2str(k(i)),'_',num2str(kmax),'.mat'];
        load(file);
        mu_list_kcat_sens_up_1000 = [mu_list_kcat_sens_up_1000;mu_list];

        file = ['kcat_sens_dn_10_',num2str(k(i)),'_',num2str(kmax),'.mat'];
        load(file);
        mu_list_kcat_sens_dn_10 = [mu_list_kcat_sens_dn_10;mu_list];

        file = ['kcat_sens_dn_1000_',num2str(k(i)),'_',num2str(kmax),'.mat'];
        load(file);
        mu_list_kcat_sens_dn_1000 = [mu_list_kcat_sens_dn_1000;mu_list];
        
        file = ['kcat_sens_up_17.47_',num2str(k(i)),'_',num2str(kmax),'.mat'];
        load(file);
        mu_list_kcat_sens_up_1747 = [mu_list_kcat_sens_up_1747;mu_list];

        file = ['kcat_sens_dn_17.47_',num2str(k(i)),'_',num2str(kmax),'.mat'];
        load(file);
        mu_list_kcat_sens_dn_1747 = [mu_list_kcat_sens_dn_1747;mu_list];
        
    end
end
cd ../../ParamSensResults;
save('mu_list_kcat_sens_up_10.mat','mu_list_kcat_sens_up_10');
save('mu_list_kcat_sens_up_1000.mat','mu_list_kcat_sens_up_1000');
save('mu_list_kcat_sens_dn_10.mat','mu_list_kcat_sens_dn_10');
save('mu_list_kcat_sens_dn_1000.mat','mu_list_kcat_sens_dn_1000');
save('mu_list_kcat_sens_up_1747.mat','mu_list_kcat_sens_up_1747');
save('mu_list_kcat_sens_dn_1747.mat','mu_list_kcat_sens_dn_1747');

end


