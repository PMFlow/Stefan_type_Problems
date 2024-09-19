%%   GRW - Stefan problem ~ Magnus
%%   ======================

clear all; close all
tic
%% Grid Initialization
I = 101; 
a=0; b=1; 
dx = (b-a)/(I-1);
x = a:dx:b;
%% Parameters
s0=dx; 
is0=2;%1; 
uD=1; %(==T_0)
D=1; % (==\alpha)
T=0.5; 
Tplot=[T/4 T/2 T]; 
past=T/3;
Dfactor=1;%2;
dtc=Dfactor*(2*D/dx^2); dt=1/dtc;
% d=1; stepU=1; U_MEAN=1; % GRW
% dt=stepU*dx/U_MEAN;
%% Initial Conditions
nT=round(T/dt); sc=zeros(1,nT);
c0=zeros(1,I); % initial temperature T(x,0)
c(1)=uD;
c(is0)=0;
q=zeros(1,I);
%% Solution
is=is0; cplot=zeros(3,I);
for it=1:nT
    t=it*dt;
    [c]=BGRW_1D(c0,I,dx,dt,q,D);
    % [c]=GRW_1D(c0,I,dx,dt,q,D,d,stepU);
    %% Boundary conditions
    c(1)=uD; % <===  Dirichlet BC left ~ [Magnus]
    c(is)=0;
    %% Diffusion front
    sc(it)=s0-(c(is)-c(is-1))/dx*dt; % (s-s0)/dt=-(c-c0)/dx ~ diffusion front ~ Ref. [21]
    is=round(sc(it)/dx);
    c0=c; s0=sc(it); c(is)=0;

    for k=1:length(Tplot)
        if  abs(t-Tplot(k))<=dt/2
            strTplot=Tplot(k);
            fprintf('t = %d\n',strTplot);
            str=['t=',num2str(strTplot)];
            strvect(k,1:length(str))=str;
            cplot(k,:)=c;
        end
    end

end
save('concentration','cplot');
%% Results
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
kt_plot=k;
plot_c(kt_plot,x,strvect); %xlim([0 2*is*dx])
iT=1:length(sc); iT=iT*dt; di=length(sc)/100;
figure
plot(iT(1:di:end),sc(1:di:end),'-+');
% semilogy(iT,sc,'-+');
xlabel('time');
ylabel('diffusion front');
t_vector=iT(1:end-1); s_GRW=sc(1:end-1);
save('sGRW','t_vector', 's_GRW');


toc
