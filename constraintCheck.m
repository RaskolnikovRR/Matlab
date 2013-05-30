%% check constraints and fval for different Q values

i = 1;
p = 0.1;
cA = [];
fvalA = [];
linA = [];

while p<1
    
    Q = [p 0 0 0 0 1-p 0 0 0 0 1-p 0 0 0 0 p]'/2;
    Ps = [0.5 0.5];
    Pxy = [p 1-p;1-p p];
    ni = [2 2 2 2];
    dels = [1 1 1 1];
    
    %% Call to convexopt.m
    
    [intmatC1,intmatC2,intmatC3] = convexopt(1,1,1,1,2,2,2,2,Pxy);
    Aeq = cat(1,intmatC1,intmatC2,intmatC3);
    beq = cat(1,transpose(Ps),[1],zeros(2*2,1));
    lin = Aeq*Q - beq;
    
    %% Call to nonlcon.m
    [c,ceq] = nonlcon(Ps,Pxy,Q,dels,ni);
    fval = costfun(Q);
    cA = cat(2,cA,c);
    fvalA = cat(2,fvalA,fval);
    linA = cat(2,linA,lin);
    p = p + 0.1;
end

display('Nonlinear constraints:');
display(cA');
display('Objective function value:');
display(fvalA');
display('Linear constraints:');
display(linA');
