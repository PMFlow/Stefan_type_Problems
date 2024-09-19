%%   diffusion in rubber 1D ~ [PEM_Surendra.pdf] - Dirichlet BC
%%   ======================
clear all; close all
tic
% iBC=1; % Dirichlet BC
iBC=2; % Robin BC
%% Grid Initialization
I = 14001; 
a=0; b=1; 
dx = (b-a)/(I-1);
x = a:dx:b;
%% Parameters
D=1e-2; beta=0.564; bb=10; H=2.5; sig=2; s0=0.01; cs0=0.5; a0=50; 
xref=10; mref=0.5; 
bb=bb/mref; cs0=cs0/mref; s0=s0/xref; 
beta=beta*xref/D; 
a0=xref*mref*a0/D; 
sig=sig*xref;
T=31*D/xref^2;
D=1;
Tplot=[T/10 T/2 T]; 
past=T/3;
Dfactor=1;
dtc=Dfactor*(2*D/dx^2); dt=1/dtc;
% d=1; stepU=1; U_MEAN=1; % GRW
% dt=stepU*dx/U_MEAN;
%% Initial Conditions
nT=round(T/dt); sc=zeros(1,nT);
c0=zeros(1,I);
q=zeros(1,I);
is0=round(s0/dx);
c0(1:is0)=cs0;
c0E=c0;c=c0;
%% Solution
is=is0; cplot=zeros(3,I);

for it=1:nT
    t=it*dt;
    [c]=BGRW_1D(c0,I,dx,dt,q,D);
    % [c]=GRW_1D(c0,I,dx,dt,q,D,d,stepU);
    %% Boundary conditions
    if iBC==1
    c(1)=uD; % <===  Dirichlet BC left(on origin)
    else
    coeff=1+beta*H*dx/D;
    c(1)=(1/coeff)*(c(2)+beta*bb*dx/D); % <===  -D*[c(2)-c(1)]/dx=\beta*(b(t)-H*c(1)) % BC left(on origin) / right(on diff. front)
    end
    c(is)=c(is-1)-(a0*dx/D)*(c(is-1)^2-c(is-1)*s0/sig); %(linearization)<=== -D*[c(is)-c(is-1)]/dx=a0*c(is)*[c(is)-\sigma];\sigma=sc/sig
    %% Diffusion front
    sc(it)=s0+a0*(c(is)-s0/sig)*dt; % diffusion front 
    is=round(sc(it)/dx);
    %
    c0=c; s0=sc(it);

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
cplot=cplot*mref;
save('concentration','cplot');
%% Results
fprintf('The space step is : %0.2e \n',dx) ;
fprintf('The time step is : %0.2e \n',dt) ;
kt_plot=k; 
plot_c(kt_plot,x,strvect); xlim([0 2*is*dx])
iT=1:length(sc); iT=iT*dt*xref^2/D/1e-2; di=2e4;
figure
plot(iT(1:di:end),xref*sc(1:di:end),'-+');
xlabel('time');
ylabel('diffusion front');

% save('I14001sig')
toc


%% main_Experiment - save('I14001sig')
% t = 3.100000e-04
% t = 1.550000e-03
% t = 3.100000e-03
% The space step is : 7.14e-05 
% The time step is : 2.55e-09 
% Elapsed time is 619.553274 seconds = 9.8614 min.

