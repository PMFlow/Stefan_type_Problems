%%   Stefan problem with kinetic condition 1D - convergence w.r.t. the analytical solution
%%   ======================
clear all; close all
tic
%% Grid Initialization
CT=0;
Itest = 1;
L2_c = zeros(4,1); eoc_c=zeros(3,1);
L2_s = zeros(4,1); eoc_s=zeros(3,1);
for iter = [100 200 400 800 1600] 
I = iter+1; 
a=0; b=1; 
dx = (b-a)/(I-1);
x = a:dx:b;
%% Parameters
alpha = 0.1; D = 1; s00 = 0.1; 
T=1; 
Tplot=[T/4 T/2 T]; 
past=T/3;
Dfactor=1;
dtc=Dfactor*(2*D/dx^2); dt=1/dtc;
%% Initial Conditions
nT=round(T/dt); sc=zeros(1,nT); sEs=zeros(1,nT); 
q=zeros(1,I);
s0=s00; mx=s0; c0=zeros(1,I); fc=zeros(1,I); 
%% manufactured solutions
Ec = @(t,x) sqrt(t)*exp(-alpha*x);
Es = @(t) s0 + 2*D*alpha*sqrt(t);
Efc = @(t,x) 1/(2*sqrt(t))*exp(-alpha*x)-D*alpha^2*sqrt(t)*exp(-alpha*x);
%% Solution
xs=0:dx:s0; is=length(xs); 
Ec0=Ec(0,xs); 
c0(1:is)=Ec0; 
for it=1:nT
    t=(it*dt)/T;
    mx=Es(t); 
    xs=0:dx:mx; Ixs=length(xs); 
    Ec0=Ec(t,xs); 
    fc(1:Ixs)=Efc(t,xs); 
    [c]=BGRW_1D(c0,I,dx,dt,q,D);
    c=c+fc*dt; 
    %% Boundary conditions
    c(1)=sqrt(t);   % left BC
    % -D*[c(is)-c(is-1)]/dx=s′(t)*c(s(t))+f2;   % right BC 
    % f2=D*alpha*c(s(t))-s′(t)*c(s(t)); ===> -D*[c(is)-c(is-1)]/dx=D*alpha*c(s(t)); ===> 
    c(is)=Ec0(is-1)-alpha*Ec0(is)*dx;
    %% Diffusion front
    % s′(t)=D*alpha(c(s(t))-alpha*s(t))+f3; 
    % f3=D*alpha/sqrt(t/T)-D*alpha(c(s(t))+alpha*s(t)); ===> 
    sc(it)=s0+D*alpha/sqrt(t)*dt;   % diffusion front
    sEs(it)=Es(t); 
    is=round(sc(it)/dx);
    %
    c0=c; s0=sc(it);

    for k=1:length(Tplot)
        if  abs(t-Tplot(k))<=dt/2
            fprintf('t = %d\n',Tplot(k));
        end
    end
    tE=t;
end
    L2_c(Itest) = ( dx )^(1/2) *norm(c(1:Ixs)-Ec0);
    L2_s(Itest) = ( dt )^(1/2) *norm(sc-sEs);
    if Itest >1
        eoc_c(Itest-1)=log10(L2_c(Itest-1)/L2_c(Itest))/log10(2);
        eoc_s(Itest-1)=log10(L2_s(Itest-1)/L2_s(Itest))/log10(2);
    end
    fprintf('The space step is : %0.2e \n',dx) ;
    fprintf('The time step is : %0.2e \n',dt) ;
    Itest = Itest+1;
toc
CT=CT+toc;
    
end
fprintf('L2_c  : %0.2e \n',L2_c)
fprintf('EOC_c : %0.2e \n',eoc_c)
fprintf('L2_s  : %0.2e \n',L2_s)
fprintf('EOC_s : %0.2e \n',eoc_s)
%% Plots
figure
ddx=5;
plot(xs(1:ddx:Ixs),c(1:ddx:Ixs),'-+',LineWidth=1); hold; 
plot(xs(1:ddx:Ixs),Ec0(1:ddx:Ixs),'-',LineWidth=1); xlim([0 is*dx]); 
legend('$c(T,x)$', '$\widetilde{c(T,x)}$','Interpreter','latex',Location='northeast',box='off');
nT=length(sc);
iT=1:nT; iT=iT*dt; di=round(nT/50);
figure
plot(iT(1:di:nT),sc(1:di:nT),'-+',LineWidth=1); hold;
plot(iT(1:di:nT),sEs(1:di:nT),'-',LineWidth=1);
xlabel('time');
ylabel('diffusion front');
legend('$s(t)$','$\widetilde{s(t)}$','Interpreter','latex',Location='northwest',box='off'); 

fprintf('total CT =  %0.2e',CT)

% main_Stef_kin_conv_new
% t = 2.500000e-01
% t = 5.000000e-01
% t = 1
% The space step is : 1.00e-02 
% The time step is : 5.00e-05 
                                    % Elapsed time is 0.329548 seconds.
% t = 2.500000e-01
% t = 5.000000e-01
% t = 1
% The space step is : 5.00e-03 
% The time step is : 1.25e-05 
                                    % Elapsed time is 1.115234 seconds.
% t = 2.500000e-01
% t = 5.000000e-01
% t = 1
% The space step is : 2.50e-03 
% The time step is : 3.13e-06 
                                    % Elapsed time is 7.127107 seconds.
% t = 2.500000e-01
% t = 5.000000e-01
% t = 1
% The space step is : 1.25e-03 
% The time step is : 7.81e-07 
                                    % Elapsed time is 44.907035 seconds.
% t = 2.500000e-01
% t = 5.000000e-01
% t = 1
% The space step is : 6.25e-04 
% The time step is : 1.95e-07 
                                    % Elapsed time is 376.947213 seconds.
% L2_c  : 1.39e-03 
% L2_c  : 4.99e-04 
% L2_c  : 1.78e-04 
% L2_c  : 6.32e-05 
% L2_c  : 2.24e-05 
                                    % EOC_c : 1.48e+00 
                                    % EOC_c : 1.49e+00 
                                    % EOC_c : 1.49e+00 
                                    % EOC_c : 1.50e+00 
% L2_s  : 1.03e-03 
% L2_s  : 5.15e-04 
% L2_s  : 2.58e-04 
% L2_s  : 1.29e-04 
% L2_s  : 6.45e-05 
                                    % EOC_s : 9.97e-01 
                                    % EOC_s : 9.98e-01 
                                    % EOC_s : 9.99e-01 
                                    % EOC_s : 1.00e+00 
% total CT =  4.30e+02
