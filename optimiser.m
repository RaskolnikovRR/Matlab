function [xA,fvalA] = optimiser()
%% Ps,Pyx definition
% Channel property Pxy and source Ps are defined
p = 0.1;
xA = [];
fvalA = [];
 
while p < 1
    Ps = [0.5 0.5];
    Pxy = [p 1-p;1-p p];
    
    %% function definition and intervals
    f = @(x,y,s,h)( (s-h)^2 + x^2);
    [xmin,ymin,hmin,smin] = deal(0);
    [xmax,smax,hmax,ymax] = deal(2);
    [xcuts,ycuts,scuts,hcuts] = deal(2);
    
    %% Call to intbymat4.m
    % to compute dels,number of elements (xi,yi,si,hi)
    [xdel,ydel,sdel,hdel,xi,yi,si,hi,iRes,F] = intbymat4(f,xmin,xmax,ymin,ymax,smin,smax,hmin,hmax,xcuts,ycuts,scuts,hcuts);
    
    %% Call to convexopt.m
    % to compute Aeq,beq,A,b
    [intmatC1,intmatC2,intmatC3] = convexopt(1,1,1,1,2,2,2,2,Pxy);
    
    %% definition of dels, ni
    dels = [xdel ydel sdel hdel];
    ni = [xi yi si hi];
    
    %% Call to fmincon.m
    Aeq = cat(1,intmatC1,intmatC2,intmatC3);
    beq = cat(1,transpose(Ps),[1],zeros(xi*yi,1));
    options = optimset('Algorithm','interior-point');
    [lb x0] = deal(zeros(xi*yi*si*hi,1));
    [x fval] = fmincon(@(Q) costfun(Q),x0,[],[],Aeq,beq,lb,[],@(Q) nonlcon(Ps,Pxy,Q,dels,ni) ,options);
    xA = cat(2,xA,x);
    fvalA = cat(1,fvalA,fval');
    p = p + 0.1;
end % main

display('Q value at minimum obj fun (each column is for a different p):');
display(xA);
display('Value of obj fun at different p:');
display(fvalA);