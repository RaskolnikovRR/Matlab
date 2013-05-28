function [c,ceq] = nonlcon(p,Q)

%........
[xdel,ydel,hdel,sdel] = deal(1);
[xi,si,hi,yi]=deal(2);
%........

[eq1,eq2,eq3,ineq2] = convexopt(p,xdel,ydel,sdel,hdel,xi,yi,si,hi,Q);
ceq = cat(1,eq1(:),transpose(eq3),eq2);
c = ineq2;
