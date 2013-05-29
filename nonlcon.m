function [c,ceq] = nonlcon(p,Q)

%% PARAMETERS, DEFINED VARIABLES
[xdel,ydel,hdel,sdel] = deal(1);
[xi,si,hi,yi]=deal(2);
p = 0.6;
Ps = [0.5 0.5];
Pxy = [p 1-p;1-p p];
%% CALL TO convexopt.m
[intmatC1,intmatC2,eq1,ineq1] = convexopt(p,xdel,ydel,sdel,hdel,xi,yi,si,hi,Q);

Aeq = cat(1,intmatC1,intmatC2);
beq = cat(1,transpose(Ps),[1]);

ceq = eq1(:);
c = ineq1;
