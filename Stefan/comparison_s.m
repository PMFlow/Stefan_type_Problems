%% compare s(t) Stefan problem

clear all; close all;

di=200;
figure; hold on
load sGRW
plot(t_vector(1:di:end), s_GRW(1:di:end),'bo-','Linewidth',1,'MarkerSize',3)
xlabel('$t$','Interpreter','latex'); 
load s1RW
plot(t_vector(1:di:end), s_RWM(1:di:end),'rs-','Linewidth',1,'MarkerSize',3)
load s1RW_ana
plot(t_vector(1:di:end), s_analytic(1:di:end),'k.-','Linewidth',1,'MarkerSize',5)
ylabel('$s(t)$','Interpreter','latex'); set(gca,'FontSize',10)
legend('GRW','RWM','Analytical','Location','southeast'); legend('Box','off')
box on;

e_RWM=norm(s_RWM-s_analytic)/norm(s_analytic)*100;
e_GRW=norm(s_GRW-s_analytic)/norm(s_analytic)*100;
fprintf('||s_RWM-s_analytic||/||s_analytic|| = %0.2e \n',e_RWM);
fprintf('||s_GRW-s_analytic||/||s_analytic||*100 = %0.2e \n',e_GRW);
