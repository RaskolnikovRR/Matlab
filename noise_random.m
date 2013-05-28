
function d = nosie_random(x,y,var)
% Input - Number of point in x and y vector in colume form: "n x 1" 
%       - var (standered variation) single value "1x1"
% Output - d :- Matrix of 10n x 2
 
n = length(x);
j = 1;
% created 10 points(noise) only for each point  
for i = 1:n
    p(j:j+9,1) = x(i,1) + var.*randn(10,1); % randome number for x
    q(j:j+9,1) = y(i,1) + var.*randn(10,1); % randome number for y
    j = j+10;
end
% Output - p for x and q for y
z=sqrt(p.^2+q.^2);
   stem3(p,q,z);
     d = [p q];
end