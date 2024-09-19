function plot_c(kt_plot,x,strvect)

load concentration

di=1; % 2; % 20; %
figure; hold all;
for k=1:kt_plot
P(k)=plot(x(1:di:end),cplot(k,1:di:end));
end
NameArray = {'Marker'}; ValueArray = {'o','+','x'}';
set(P,NameArray,ValueArray);
xlabel('$x$','Interpreter','latex');
ylabel('$c(x,t)$','Interpreter','latex');
legend(strvect); legend('boxoff'); % xlim([0 100])

