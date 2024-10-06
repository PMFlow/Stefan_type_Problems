function [fs]=Fs(t,T,cs,s0,a0,sig)

    % % with Matlab:
    % syms s(t);
    % syms cs;
    % syms s0;
    % syms a0;
    % syms sig;
    % syms T;
    % s = t/T/10-t^2/T^2/10+t^3/T^3/30+s0;
    % fs = diff(s,t) - a0*(cs-s/sig);
    
    fs = s0*(1/(10*T) - t/(5*T^2) + t^2/(10*T^3)) - a0*(cs - (s0*(t/(10*T) - t^2/(10*T^2) + t^3/(30*T^3) + 1))/sig);
