%triangle pulse function
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: pulse evaluated at t
function res = triangle_pulse(t,w,h)
    t = t*(2/w);
    res = 1-min(1*abs(t-1),1);
    res = h*res;
end
