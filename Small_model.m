
mu_list = 0.1:0.025:0.625;

flux = [];
rel_flux = [];

for i = 1:length(mu_list)
    mu = mu_list(i);

    A = [-3 -2 -2;0.0073 0.002 0.0054];
    
    b = [-80.44*mu+1.64 0.05];

    Aeq = [];
    beq = [];

    lb = [0 0 0];
    ub = [100 100 1.6568*mu];

    f = [1 1 0];

    x = linprog(f,A,b,Aeq,beq,lb,ub);
    
    flux = [flux x];
    rel_flux = [rel_flux [x(1)/(x(1)+x(2));x(2)/(x(1)+x(2))]];
    
end

clear i mu A b Aeq beq lb ub f x factor;

figure('Name','1');
hold on;
plot(mu_list,rel_flux(1,:),'-o','LineWidth',1.5,'Color',[55,126,184]/256,'MarkerSize',8,'MarkerEdgeColor',[55,126,184]/256,'MarkerFaceColor',[55,126,184]/256);
plot(mu_list,rel_flux(2,:),'-o','LineWidth',1.5,'Color',[228,26,28]/256,'MarkerSize',8,'MarkerEdgeColor',[228,26,28]/256,'MarkerFaceColor',[228,26,28]/256);
ylabel('Fraction in total glucose carbon','FontSize',12,'FontName','Helvetica');

set(gca,'FontSize',10,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',12,'FontName','Helvetica');
legend({'Mixed acid','Lactate'},'FontSize',12,'FontName','Helvetica','location','west');
set(gcf,'position',[0 0 335 300]);
set(gca,'position',[0.11 0.11 0.75 0.8]);

figure('Name','2');
plot(mu_list,flux(3,:),'-o','LineWidth',1.5,'Color',[77,175,74]/256,'MarkerSize',8,'MarkerEdgeColor',[77,175,74]/256,'MarkerFaceColor',[77,175,74]/256);
set(gca,'FontSize',10,'FontName','Helvetica');
ylabel('q(mmol/gCDW/h)','FontSize',12,'FontName','Helvetica');
xlabel('Growth rate (/h)','FontSize',12,'FontName','Helvetica');
legend({'Arginine'},'FontSize',12,'FontName','Helvetica','location','nw');
set(gcf,'position',[0 0 335 300]);
set(gca,'position',[0.12 0.11 0.75 0.8]);