function plot_c_conv(kt_plot,strvect)

load concentration_conv

figure; hold all;
for k=1:kt_plot
P(k)=plot(xk(k,:),cplot(k,:));
end
NameArray = {'Marker'}; ValueArray = {'o','+','x'}';
set(P,NameArray,ValueArray);
xlabel('$x$','Interpreter','latex');
ylabel('$c(x,t)$','Interpreter','latex');
legend(strvect); legend('boxoff'); 

