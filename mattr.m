function [ress] = mattr(x,v)

% f(x,v) = @(x,v) (exp(-x^2/(2*s^2))*exp(-v^2/2)/(2*pi*s)*(k^2*(s*signum(x) - x)^2 + ( s*signum(x) - s*tanh( s*(s*signum(x) + v)))^2 ) );
s = 1;
k = 1;

% ress = f(x,v);
syms i j ;
i = 1;
j = 1;
% ress= ones(length(x),length(v));
while i <= length(x)
    j = 1;
    while j <= length(v)
        ress(i,j) = vpa(exp(-x(i)^2/(2*s^2))*exp(-v(j)^2/2)/(2*pi*s)*(k^2*(s*signum(x(i)) - x(i))^2 + ( s*signum(x(i)) - s*tanh( s*(s*signum(x(i)) + v(j))))^2 ) );
        j = j+1;
    end
    i = i+1;
end
ress