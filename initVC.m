function[xx]=initVC(m,k,u,v)

syms t s X;

X=(m*s*u+m*v)/(m*s^2+k); %this is a simpler code because taking laplace and solving steps have been done by user; HERE xx=X(s)

xx=ilaplace(X); %take inverse to get x(t)
