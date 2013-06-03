function [intmats,intmatC1,intmatC2,intmatC3] = convexopt(xdel,ydel,sdel,hdel,xi,yi,si,hi,Pxy)
%% PARAMETERS AND VARIABLE DEFINITIONS
syms rowiter li;
li = xi*yi*si*hi;

%% CONSTRAINT1 : INTEGRATION OF Q over x,y,h = Ps(S)
% int(x,y,h)Q = Ps(s)
rowiter = 1;
intmat1 = zeros(li/xi,li);

while( rowiter <= li/xi )
    intmat1(rowiter,(xi*(rowiter-1)+1):xi*rowiter) =  1;
    rowiter = rowiter + 1;
end

% intmat1 = (li/xi,li)

rowiter = 1;
intmat2 = zeros(li/(xi*yi),li/xi);

while(rowiter <= li/(xi*yi))
    intmat2(rowiter,(yi*(rowiter-1)+1):yi*rowiter) = 1;
    rowiter = rowiter + 1;
end

intmat3 = zeros(si,si*hi);
rowiter = 1;
while rowiter <= si
    coliter = 1;
    while( coliter <= hi)
        intmat3(rowiter,rowiter + si*(coliter - 1)) = 1;
        coliter = coliter + 1;
    end
    rowiter = rowiter + 1;
end

intmatC1 = (xdel*ydel*hdel)*intmat3*intmat2*intmat1;
% intmatC1 (si*li) = (si,si*hi)*(si*hi,si*hi*yi)*(si*hi*yi.li)
% Aeq*Q = Ps = beq

%% CONSTRAINT2 - (x,y,h,s) Q  = 1

intmatC2 = xdel*ydel*hdel*sdel*ones(1,li);
% intmatC2 (1,li)
% Aeq*Q = beq = 1

%% CONSTRAINT3 : INTEGRATION OF Q over x,y = P(Y|X)*b(X)
% int(s,h) Q  = P(y/x)*b(x) for all y,x
% P(Y|X) is a matrix: (xi,yi)
% step1 : rearrange Q in order X -> Y -> S -> H

rowiter = 1;
intmat4 = zeros(li/hi,li);
while(rowiter <= li/hi)
    coliter =1;
    while(coliter <= hi)
        intmat4(rowiter,rowiter + (coliter - 1)*li/hi) = 1;
        coliter = coliter + 1;
    end
    rowiter = rowiter + 1;
end

rowiter = 1;
intmat5 = zeros(li/hi/si,li/hi);
while(rowiter <= li/hi/si)
    coliter = 1;
    while( coliter <= si)
        intmat5(rowiter,rowiter + (coliter - 1)*li/hi/si) = 1;
        coliter = coliter + 1;
    end
    rowiter = rowiter + 1;
end

rowiter = 1;
intmat6 = zeros(li/hi/si/yi,li/hi/si);
while(rowiter <= li/hi/si/yi)
    coliter = 1;
    while(coliter <= yi)
        intmat6(rowiter,rowiter + (coliter - 1)*li/hi/si/yi) = 1;
        coliter = coliter + 1;
    end
    rowiter = rowiter + 1;
end

intmatC3a = (ydel*sdel*hdel)*(intmat6*intmat5*intmat4);
%( li/hi/si/yi,li)  = (li/hi/si/yi,li/hi/si)*(li/hi/si,li/hi) * (li/hi,li)
% (xi,li)

vecPxy = Pxy(:);
diagPxy = zeros(xi*yi,xi*yi);
diagPxy = diag(vecPxy,0);

repBX = repmat( diag(ones(1,xi),0),yi,1);
intmatC3 = diagPxy*repBX*intmatC3a;
% (xiyi,li) = (xiyi,xiyi)*(xiyi,xi)*(xi,li)

intmatC3 = intmatC3 - intmat5*intmat4;

%% pass back intmats
intmats = {intmat1 intmat2 intmat3 intmat4 intmat5 intmat6};

%% LOCAL FUNCTIONS

    function rQ = rearrange(Q,xi,yi,si,hi)
        
        syms X Y S H i; % counters
        
        X = 1;
        i = 1;
        while (X <= xi)
            Y = 1;
            while(Y <= yi)
                H = 1;
                while(H <= hi)
                    S = 1;
                    while(S <= si)
                        now = xi*yi*(S-1) + (H-1)*xi*yi*si + xi*(Y-1) + X;
                        rQ(i) = Q(now);
                        S = S + 1;
                        i = i + 1;
                    end
                    H = H + 1;
                end
                Y = Y + 1;
            end
            X = X + 1;
        end
        
    end

    function v = vecinv(vec)
        
        i = 1;
        while( i <= length(vec))
            v(i,:) = 1/vec(i)/length(vec);
            i = i + 1;
        end
        
        
    end
end % main function

