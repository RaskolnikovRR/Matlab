function[num,den] = transferFunc(m,k,u,v)
syms m k u v X;
syseq='m*D(D(x))(t)+k*x(t)';
lapSys=laplace(syseq,t,s);
lapSys=subs(lapSys,{'laplace(x(t),t,s)','D(x)(0)','x(0)',{X,v,u});
transF=solve(lapSys,X);
finalX=ilaplace(transF);
finalX
