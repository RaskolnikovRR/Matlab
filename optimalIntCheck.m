function [iResult] = optimalIntCheck(xlm,ylm,s,k)
timestamp
syms x v;
func = @(x,v) (exp(-x.^2/2/s^2).*exp(-v.^2/2)/2/pi/s).*( k^2*(s*signum(x) - x).^2 + ( s*signum(x) - s*tanh(s*(s*signum(x) + v)) ).^2 );
iResult = dblquad(func,xlm*-1,xlm,ylm*-1,ylm);
timestamp
