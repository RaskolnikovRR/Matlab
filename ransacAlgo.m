function [Etotal,a,b,Ea,Eb] = ransacAlgo(n, iterations)


format short;
syms data;

% generate data points
data = rand(2,n);


syms x y a b;
ErrorEq = ( y - a*x + b ); %/ sqrt(a^2 + 1);

% start main iteration for a,b
syms Etotal Ea Eb;

%for j = 1:iterations
    
    Etotal = 0;
    
    for i = 1:n
        Ei = subs(ErrorEq, {y,x},{data(2,i),data(1,i)});
        Etotal = Etotal + Ei/n;
    end
    
    Etotal = simplify(Etotal);
    Etotal = vpa(Etotal , 5);
    Ea = diff(Etotal,'a');
    Ea = vpa(Ea,5);
    Eb = diff(Etotal,'b');
    Eb = vpa(Eb,5);
    Eq1 = strcat(char(Ea),'=0');
    Eq2 = strcat(char(Eb),'=0');
    [a,b] = solve(Eq1,Eq2);
%     a = a( real(a) == a);
%     b = b( real(b) == b);
     
    
%end


   
    