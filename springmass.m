function [response]= springmass(m,k,b,f)

syms x t s Xs;
%writing the system equation from fbd
syseq='f-k*x(t)-b*D(x)(t)=m*D(D(x))(t)'

%substituting values of m,k,b,f
syseq=subs(syseq,{'m','k','b','f'},{m,k,b,f})

%laplace of the system equation
lapSys=laplace(syseq,t,s)
lapSys=subs(lapSys,{'laplace(x(t),t,s)','D(x)(0)','x(0)'},{Xs,0,0}) %intial conditions zero
Xs=solve(lapSys,Xs);
response=ilaplace(Xs);