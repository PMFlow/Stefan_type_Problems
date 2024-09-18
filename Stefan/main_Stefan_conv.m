%%   GRW - Stefan problem ~ Magnus - grid-convergence test
%%   ======================

clear all; close all
%% Grid Initialization
CT=0;
Itest = 1;
L2_c = zeros(5,1); eoc_c=zeros(4,1);
L2_s = zeros(5,1); eoc_s=zeros(4,1);
for iter = [25 50 100 200 400] % [50 100 200 400 800] % [100 200 400 800 1600] % 
tic
    I = iter+1;
    a=0; b=1;
    dx = (b-a)/(I-1);
    x = a:dx:b;
    %% Parameters
    s0=dx;
    is0=2;%1;
    alpha=1; % K/(rho*c); % Thermal diffusivity.
    beta=1; % l/c; % Parameter with unit [K].
    T_0=1; % [degree C] Temperature for the constant temperature BC.
    T=0.5;
    Tplot=[T/4 T/2 T];
    past=T/3;
    Dfactor=1; % with BGRW
    dtc=Dfactor*(2*alpha/dx^2); dt=1/dtc;
    % d=1; stepU=1; U_MEAN=1; % with unbiased GRW
    % dt=stepU*dx/U_MEAN;
    %% Initial Conditions
    nT=round(T/dt); sc=zeros(1,nT);
    c0=zeros(1,I); % initial temperature T(x,0)
    c(1)=T_0;
    c(is0)=0;
    q=zeros(1,I);
    %% Solution
    is=is0; cplot=zeros(3,I);
    for it=1:nT
        t=it*dt;
        [c]=BGRW_1D(c0,I,dx,dt,q,alpha);
        % [c]=GRW_1D(c0,I,dx,dt,q,alpha,d,stepU);
        %% Boundary conditions
        c(1)=T_0; % <===  Dirichlet BC left ~ Ref. [21]
        c(is)=0;
        %% Diffusion front
        sc(it)=s0-(c(is)-c(is-1))/dx*dt; % (s-s0)/dt=-(c-c0)/dx ~ diffusion front ~ Ref. [21]
        is=round(sc(it)/dx);
        c0=c; s0=sc(it); c(is)=0;

        for k=1:length(Tplot)
            if  abs(t-Tplot(k))<=dt/2
                strTplot=Tplot(k);
                str=['t=',num2str(strTplot)];
                strvect(k,1:length(str))=str;
                cplot(k,:)=c;
            end
        end

    end
    %% Results
    s_GRW=sc(1:end-1); c_GRW=c;
    iT=1:length(sc); iT=iT*dt; 
    % For constant temperature or special heat flux BC, the analytical solution requires the solution of the transcendental equation:
    lambda=trans_eq(beta,T_0); 
    c_ana_vector = T_0*(1-(erf(x/(2*sqrt(alpha*t))))/(erf(lambda)));
    for i=1:I
        if c_ana_vector(i) < 0
            c_ana_vector(i)=0;
        end
    end
    s_ana_vector = dx+2*lambda*sqrt(alpha*iT(1:end-1));
    L2_c(Itest) = ( dx )^(1/2) *norm(c_GRW-c_ana_vector);
    L2_s(Itest) = ( dt )^(1/2) *norm(s_GRW-s_ana_vector);
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
save('concentration','cplot');
kt_plot=k;
plot_c(kt_plot,x,strvect); 
di=length(sc)/100;
figure
plot(iT(1:di:end),sc(1:di:end),'-+');

fprintf('total CT =  %0.2e',CT)

function [x0] = trans_eq(beta,T_0)
f = @(x) sqrt(pi)*beta*x*exp(x^2)*erf(x)-T_0;
fprim = @(x) beta*(sqrt(pi)*exp(x^2)*erf(x)*(2*x^2+1)+2*x^2);
tol=1e-6;
x0=1;
while abs(f(x0)) > tol
    x0 = x0-f(x0)/fprim(x0);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% convergence of c & s: 4 successive halvings of dx starting with dx=4e-2; 
% main_Stefan_conv
% The space step is : 4.00e-02 
% The time step is : 8.00e-04 
% Elapsed time is 0.076764 seconds.
% The space step is : 2.00e-02 
% The time step is : 2.00e-04 
% Elapsed time is 0.066054 seconds.
% The space step is : 1.00e-02 
% The time step is : 5.00e-05 
% Elapsed time is 0.127952 seconds.
% The space step is : 5.00e-03 
% The time step is : 1.25e-05 
% Elapsed time is 0.532853 seconds.
% The space step is : 2.50e-03 
% The time step is : 3.13e-06 
% Elapsed time is 3.011082 seconds.
% L2_c  : 6.52e-03 
% L2_c  : 3.49e-03 
% L2_c  : 1.84e-03 
% L2_c  : 9.22e-04 
% L2_c  : 4.28e-04 
                                % EOC_c : 9.00e-01 
                                % EOC_c : 9.24e-01 
                                % EOC_c : 9.98e-01 
                                % EOC_c : 1.11e+00 
% L2_s  : 9.62e-03 
% L2_s  : 4.89e-03 
% L2_s  : 2.41e-03 
% L2_s  : 1.18e-03 
% L2_s  : 5.76e-04 
                                % EOC_s : 9.77e-01 
                                % EOC_s : 1.02e+00 
                                % EOC_s : 1.03e+00 
                                % EOC_s : 1.03e+00 
% total CT =  3.82e+00>> 

