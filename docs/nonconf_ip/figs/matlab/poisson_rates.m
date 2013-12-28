set(0,'defaultlinelinewidth',1.5)

h =  [0.25, 0.125, 0.0625, 0.03125]
conf_err =  [0.11292510250062097, 0.031661312420476179, 0.0081570955718738691, 0.0020549140931555316]
nconf_err =  [0.049155516397130262, 0.012556193865282748, 0.0031450644387712398, 0.00078564248820516631]

loglog(h,conf_err,'b-');hold on;loglog(h,nconf_err,'bo-')

h =  [0.25, 0.125, 0.0625, 0.03125]
conf_err =  [0.0047368635504671094, 0.00056317792376020276, 6.9233636338376692e-05, 8.616169151022162e-06]
nconf_err =  [0.0034947390312099623, 0.00039451618850192769, 4.7546196697389755e-05, 5.8744142074609405e-06]

loglog(h,conf_err,'r-');hold on;loglog(h,nconf_err,'ro-')

conf_err =  [0.00035302674748083845, 2.1792879237613943e-05, 1.3431462980237558e-06, 8.32272311361012e-08]
nconf_err =  [0.00033725807257816418, 2.1122627031359447e-05, 1.3092841843400464e-06, 8.1329374119981046e-08]

loglog(h,conf_err,'k-');hold on;loglog(h,nconf_err,'ko-')

legend('H^1 p = 1','Broken p = 1',...
    'H^1 p = 2','Broken p = 2',...
    'H^1 p = 3','Broken p = 3')
xlabel('Mesh size h','fontsize',14)
ylabel('L^2 error','fontsize',14)
set(gca,'fontsize',14)
%set(gcf,'DefaultAxesFontSize', 12)
%set(gcf,'DefaultTextFontSize', 12)
%axes('FontSize',14)

%print(gcf,'-depsc','../poisson_rates.eps')
%print(gcf,'-dpng','../poisson_rates.png')
