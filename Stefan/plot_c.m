function plot_c(kt_plot,x,strvect)

load concentration

di=1; % 10; % 20;
figure; hold all;
for k=1:kt_plot
P(k)=plot(x(1:di:end),cplot(k,1:di:end),'MarkerSize',3);
end
NameArray = {'Marker'}; ValueArray = {'o','+','x'}';
set(P,NameArray,ValueArray);
xlabel('$x$','Interpreter','latex');
ylabel('$m(x,t)$','Interpreter','latex');
legend(strvect); legend('boxoff'); xlim([0 1.]); box on;

