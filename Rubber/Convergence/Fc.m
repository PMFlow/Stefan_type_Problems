function [fc,Eq2,Eq3]=Fc(t,x,T,beta,b,H,cs,s0,a0,sig)

    %% with Matlab:
    % syms c(t,x);
    % syms mx;
    % syms T;
    % syms beta;
    % syms b;
    % syms H;
    % syms cs;
    % syms s0;
    % syms a0;
    % syms sig;
    % s = s0*(t/T/10-t^2/T^2/10+t^3/T^3/30+1);    
    % c=(1-x/s).^3*cos(t/T);
    % fc = diff(c,t)-diff(diff(c,x),x);  
    % Eq2 = -diff(c,x)-beta*(b-H*c);
    % Eq3 = -diff(c,x)-a0*(cs-s/sig)*c; 

    fc = (sin(t/T)*(x/(s0*(t/(10*T) - t^2/(10*T^2) + t^3/(30*T^3) + 1)) - 1).^3)/T + (6*cos(t/T)*(x/(s0*(t/(10*T) - t^2/(10*T^2) ...
         + t^3/(30*T^3) + 1)) - 1))/(s0^2*(t/(10*T) - t^2/(10*T^2) + t^3/(30*T^3) + 1)^2) + (3*x*cos(t/T).*(x/(s0*(t/(10*T) - t^2/(10*T^2) ...
         + t^3/(30*T^3) + 1)) - 1).^2*(1/(10*T) - t/(5*T^2) + t^2/(10*T^3)))/(s0*(t/(10*T) - t^2/(10*T^2) + t^3/(30*T^3) + 1)^2);

    Eq2 = 3*cos(t/T)/(s0*(t/(10*T) - t^2/(10*T^2) + t^3/(30*T^3) + 1)) - beta*b + beta*H*cos(t/T);
    Eq3 = 0.*x;







