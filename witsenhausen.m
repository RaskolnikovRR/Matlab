function [J_bar] = witsenhausen(f1,f2,sigma,k)

% f1 f2 ==> user input of functions f,g 
% example: for f = sin, input f1 should be '@(x) sin(x)'

% witsenhausen cost = J(f,g,x,v) =
% k^2*(f(x) - x)^2 + (f(x) - g( f(x) + v) )^2
% x = N(0,sigma^2)
% v = N(0,1)
% x,v are independent random variables
% J is a function of the two random variables X,V which are INDEPENDENT

% E( J(X,V) ) = intX(-inf --> inf) intY( -inf --> inf) J(x,v) f(x,v) dx dv
% Since x,v are independent joint probability distribution function -
% - f(x,v) => fx(X)fv(V)

syms x v J
% gaussian distributions x,v
fx = @(x) exp( -x^2/2/sigma^2)/sqrt(2*pi*sigma^2);
fv = @(v) exp( -v^2/2)/sqrt(2*pi);

% constructing functions f,g from strings f1,f2
f = str2func(f1);
g = str2func(f2);

% defining cost function J
J = @(x,v) (k^2*(f(x) - x)^2 + (f(x) - g(f(x) + v) )^2)';

% integration
toIntegrate = vpa(J(x,v)*fx(x)*fv(v),2);
int1 = vpa(int(toIntegrate,v,-inf,inf),2);
J_bar = vpa(int(int1,x,-inf,inf),2);
