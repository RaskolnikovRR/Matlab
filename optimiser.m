function [xA,fvalA] = optimiser()
%% function definition and intervals
f = @(x,y,s,h)( (s-h)^2 + x^2);
[xmin,ymin,hmin,smin] = deal(0);
[xmax,smax,hmax,ymax] = deal(12);
[xcuts,ycuts,scuts,hcuts] = deal(12);

%% del calculation
xdel = (xmax - xmin)/xcuts;
ydel = (ymax - ymin)/ycuts;
sdel = (smax - smin)/scuts;
hdel = (hmax - hmin)/hcuts;

%% Ps,Pyx definition
% Channel property Pxy and source Ps are defined
sig = 1;
Ps = [smin + sdel:sdel:smax];
Ps = exp(-Ps.^2/2/sig^2)/sig/sqrt(2*pi);

Pxy = zeros(xcuts,ycuts);
curX = xmin;
curY = ymin;

ii = 1;
check = 0;
while ii <= ycuts
    jj = 1;
    curX = xmin;
    while jj <= xcuts
        Pxy(jj,ii) = vpa(exp( -(curY - curX)^2/2)/sqrt(2*pi),32);
        curX = curX + xdel;
        jj = jj + 1;
        check = check + 1;
    end
    curY = curY + ydel;
    ii = ii + 1;
end

%% Call to intbymat4.m
% to compute dels,number of elements (xi,yi,si,hi)
[xdel,ydel,sdel,hdel,xi,yi,si,hi,iRes,F] = intbymat4(f,xmin,xmax,ymin,ymax,smin,smax,hmin,hmax,xcuts,ycuts,scuts,hcuts);

%% Call to convexopt.m
% to compute Aeq,beq,A,b
[intmats,intmatC1,intmatC2,intmatC3] = convexopt(xdel,ydel,sdel,hdel,xi,yi,si,hi,Pxy);

%% definition of dels, ni
dels = [xdel ydel sdel hdel];
ni = [xi yi si hi];

%% Call to fmincon.m
Aeq = cat(1,intmatC1,intmatC2,intmatC3);
beq = cat(1,transpose(Ps),[1],zeros(xi*yi,1));
options = optimset('Algorithm','interior-point');
[lb x0] = deal(zeros(xi*yi*si*hi,1));
[x fval] = fmincon(@(Q) costfun(Q),x0,[],[],Aeq,beq,lb,[],@(Q) nonlcon(Ps,Pxy,Q,dels,ni,intmats) ,options);
xA = cat(2,xA,x);
fvalA = cat(1,fvalA,fval');
nonlcon(Ps,Pxy,x,dels,ni,intmats)

display('Q value at minimum obj fun (each column is for a different p):');
display(xA);
display('Value of obj fun at different p:');
display(fvalA);