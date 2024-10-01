%%   diffusion in rubber 1D - convergence w.r.t. the manufactured solution
%%   ======================
clear all; close all
tic
% iBC=1; % Dirichlet BC 
iBC=2; % Robin BC
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
beta=1; bb=1; s00=0.1; H=1; sig=1; a0=1;
T=1e-3; 
Tplot=[T/4 T/2 T]; 
past=T/3;
D=1;
Dfactor=1;
dtc=Dfactor*(2*D/dx^2); dt=1/dtc;
%% Initial Conditions
nT=round(T/dt); sc=zeros(1,nT); sEs=zeros(1,nT); 
q=zeros(1,I);
s0=s00; mx=s0; c0=zeros(1,I); fc=zeros(1,I); 
%% manufactured solutions
Ec = @(t,x) (1-x/mx).^3*cos(t/T);
Es = @(t) t/T/10-t^2/T^2/10+t^3/T^3/30+s00;
%% Solution
xs=0:dx:s0; is=length(xs); 
Ec0=(1-xs/mx).^3; 
c0(1:is)=Ec0; 
for it=1:nT
    t=it*dt;
    mx=Es(t); 
    xs=0:dx:mx; Ixs=length(xs);
    Ec0=(1-xs/mx).^3*cos(t/T); 
    [Efc,Eq2,Eq3]=Fc(t,xs,T,beta,b,H,0,s00,a0,sig); 
    fc(1:Ixs)=Efc;    
    [fs]=Fs(t,T,0,s00,a0,sig);  
    [c]=BGRW_1D(c0,I,dx,dt,q,D);
    c=c+fc*dt; 
    %% Boundary conditions
    if iBC==1
        c(1)=Ec0(1); % Dirichlet BC left(on origin) 
    else
        % -D*[c(2)-c(1)]/dx=\beta*(b(t)-H*c(1)) % Robin BC left(on origin) ===> 
        coeff=1+beta*H*dx/D;
        c(1)=(1/coeff)*(Ec0(2)+beta*bb*dx/D)+Eq2(1)*dx/D; 
    end
    % -D*[c(is)-c(is-1)]/dx=a0*c(is)*[c(is)-\sigma]+fs;  \sigma=sc/sig (explicit linearization, is ---> is-1) ===> 
    c(is)=Ec0(is-1)-(a0*dx/D)*(Ec0(is)^2-Ec0(is)*s0/sig)-Eq3(is)*dx/D; 
    %% Diffusion front
    sc(it)=s0+a0*(c(is)-s0/sig)*dt+fs*dt; % diffusion front
    sEs(it)=Es(t); 
    is=round(sc(it)/dx);
    %
    c0=c; s0=sc(it);

    for k=1:length(Tplot)
        if  abs(t-Tplot(k))<=dt/2
            strTplot=Tplot(k);
            fprintf('t = %d\n',strTplot);
            str=['t=',num2str(strTplot)];
            strvect(k,1:length(str))=str;
            xk(k,1:Ixs)=xs;
            cplot(k,1:Ixs)=c(1:Ixs);
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
save('concentration_conv','xk','cplot');
kt_plot=k;
plot_c_conv(kt_plot,strvect); hold all;
plot(xs,Ec0,'-',LineWidth=1); xlim([0 is*dx])
nT=length(sc);
iT=1:nT; iT=iT*dt; di=round(nT/50);
figure
plot(iT(1:di:nT),sc(1:di:nT),'-+'); hold;
plot(iT(1:di:nT),sEs(1:di:nT),'-',LineWidth=1);
% semilogy(iT,sc,'-+');
xlabel('time');
ylabel('diffusion front');
legend('$s(t)$','$\widetilde{s(t)}$','Interpreter','latex',Location='northwest',box='off'); 

fprintf('total CT =  %0.2e',CT)

%% main_RobinBC_conv.m: Robin BC; 4 successive halvings starting with dx=1e-2; 
% main_RobinBC_conv
% t = 2.500000e-04
% t = 5.000000e-04
% t = 1.000000e-03
% The space step is : 1.00e-02 
% The time step is : 5.00e-05 
                                    % Elapsed time is 0.045562 seconds.
% t = 2.500000e-04
% t = 5.000000e-04
% t = 1.000000e-03
% The space step is : 5.00e-03 
% The time step is : 1.25e-05 
                                    % Elapsed time is 0.057908 seconds.
% t = 2.500000e-04
% t = 5.000000e-04
% t = 1.000000e-03
% The space step is : 2.50e-03 
% The time step is : 3.13e-06 
                                    % Elapsed time is 0.083846 seconds.
% t = 2.500000e-04
% t = 5.000000e-04
% t = 1.000000e-03
% The space step is : 1.25e-03 
% The time step is : 7.81e-07 
                                    % Elapsed time is 0.229305 seconds.
% t = 2.500000e-04
% t = 5.000000e-04
% t = 1.000000e-03
% The space step is : 6.25e-04 
% The time step is : 1.95e-07 
                                    % Elapsed time is 0.936351 seconds.
% L2_c  : 1.95e-03 
% L2_c  : 4.33e-04 
% L2_c  : 1.04e-04 
% L2_c  : 2.54e-05 
% L2_c  : 6.31e-06 
                                    % EOC_c : 2.17e+00 
                                    % EOC_c : 2.06e+00 
                                    % EOC_c : 2.03e+00 
                                    % EOC_c : 2.01e+00 
% L2_s  : 5.81e-05 
% L2_s  : 1.44e-05 
% L2_s  : 3.61e-06 
% L2_s  : 9.01e-07 
% L2_s  : 2.25e-07 
                                    % EOC_s : 2.01e+00 
                                    % EOC_s : 2.00e+00 
                                    % EOC_s : 2.00e+00 
                                    % EOC_s : 2.00e+00 
% total CT =  1.35e+00
