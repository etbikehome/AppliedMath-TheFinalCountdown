%INPUTS:
%x: vector of points to evaluate mode shape
%n: index of mode to calculate
%L: length of string
%c: wave speed
%OUTPUTS:
%V: mode shape y values for every element of x
%freq: resonant frequency corresponding to this mode
function [V, freq] = modes(x, n, L, c)
    V = sin(pi*n/L*x);
    freq = c*pi*n/L;
end