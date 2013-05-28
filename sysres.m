% system response for different inputs
function [ c ]=sysres(r)

syms s g t;
% c,r in frequency domain
g = [2 1 3 1 3]*[s^4 ;s^3; s^2; s; 1];

% g = transfer function
c = g*r;
CinTIME = ilaplace(c)
ezplot(CinTIME);