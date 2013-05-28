function[finalX] = initV(m,k,u,v)

syms X t s;
syseq='m*D(D(x))(t)+k*x(t)'; %starts from the system equation

syseq=subs(syseq,{'m','k'},{m,k}); %substitute m,k values

lapSys=laplace(syseq,t,s) %taking laplace of entire equation

lapSys=subs(lapSys,{'laplace(x(t),t,s)','D(x)(0)','x(0)'},{X,v,u}) %substitute X

XF=solve(lapSys,X); %solve equation to get X(s)

finalX=ilaplace(XF); %take inverse to get x(t)
