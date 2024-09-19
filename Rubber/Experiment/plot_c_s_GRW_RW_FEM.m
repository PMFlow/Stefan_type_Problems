%% plots comparison of c(t) and s(t) solutions obtained by different methods
clear all; close all

load I14001sig.mat 

di=100; %
figure; hold on; % concentration
plot(xref*x(1:di:end),cplot(3,1:di:end),'b-o','LineWidth',1,'MarkerSize',3); xlim([0 2.5]);
load matlab_data_con_RW_rubber1.mat %% load random walk solution (concentration profile at T = 31 minutes)
plot(xx(25:5:end), c_RW(25:5:end),'rs-','LineWidth',0.75,'MarkerSize',3);
load matlab_data_con_FEM_rubber1.mat % load finite element solution (concentration profile at T = 31 minutes)
plot(xx(1:5:end), c_t_x(1:5:end),'k.-','LineWidth',1, 'MarkerSize',5);
xscale log; set(gca, 'XTick', [0.1 0.5 1 2.5]);
xlabel('$x$','Interpreter','latex');
ylabel('$c(x,t)$','Interpreter','latex'); box on;
legend('GRW','RWM','FEM','Location','northeast','box','off'); 

iT=1:length(sc); iT=iT*dt*xref^2/D/1e-2;
di=20000; %
figure; hold on; % diffusion front
plot(iT(1:di:end),xref*sc(1:di:end),'b-o','LineWidth',1,'MarkerSize',3);
load matlab_data_st_RW_rubber1.mat % load random walk solution (moving boundary)
plot(t(1:1000:end),s_t(1:1000:end),'r-s','LineWidth',0.75,'MarkerSize',3);
load matlab_data_st_FEM_rubber1.mat % load finite element solution (moving boundary)
plot(t(1:di:end),s_t(1:di:end),'k.-','LineWidth',1,'MarkerSize',5);
time = [0, 3.5, 10, 30];
s =  [0, 1, 2, 2] ;
plot(time,s,'k*','MarkerSize',5);
xlabel('time');
ylabel('diffusion front'); box on;
legend('GRW','RWM','FEM','Experiemnt','Location','southeast','box','off'); 
