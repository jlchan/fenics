% dp 1
h =  [0.125]
energy_err =  [0.22858090993977098];
L2_err =  [0.012003627548950948];

% dp 2
energy_err = [energy_err 0.23029652417346533];
L2_err =  [L2_err 0.011963812622815932];

% dp 3
energy_err =  [energy_err 0.23078143548994856];
L2_err =  [L2_err 0.011953146208037663];

% dp 4
energy_err = [energy_err 0.23090056223942759]
L2_err = [L2_err .011950356327737074]

figure
plot(1:length(L2_err), energy_err,'bs-')
hold on
plot(1:length(L2_err), L2_err,'ro-')
set(0,'defaultlinelinewidth',1.5)
xlabel('dp','fontsize',15)
legend('L^2 error', 'Energy error')
set(gca,'fontsize',14)
% print(gcf,'-depsc','../dp_eps1e0.eps')


% dp 1
energy_err =  [0.13040670636936783];
L2_err =  [0.12676314256288002];

% dp 2
energy_err = [energy_err 0.1304088510666537];
L2_err = [L2_err 0.12676333852409472];

% dp 3
energy_err = [energy_err 0.13041019002765655];
L2_err = [L2_err 0.12676338950353838];

% dp 4
energy_err = [energy_err 0.13041096408992561]
L2_err =  [L2_err 0.12676342243657773]

figure
set(0,'defaultlinelinewidth',1.5)
plot(1:length(L2_err), energy_err,'bs-')
hold on
semilogy(1:length(L2_err), L2_err,'ro-')
xlabel('dp','fontsize',15)
set(gca,'fontsize',14)
legend('L^2 error', 'Energy error')
% print(gcf,'-depsc','../dp_eps1e4.eps')