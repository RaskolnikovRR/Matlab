function [x fval] = witsenhausen(small,large,gap)

tic;
%% function definition and intervals
[xmin,ymin,hmin,smin] = deal(small);
[xmax,smax,hmax,ymax] = deal(large);
[xcuts,ycuts,scuts,hcuts] = deal(gap);

cutscell = num2cell([xcuts ycuts scuts hcuts]);
[xi yi si hi] = cutscell{:};

%% del calculation
xdel = (xmax - xmin)/xcuts;
ydel = (ymax - ymin)/ycuts;
sdel = (smax - smin)/scuts;
hdel = (hmax - hmin)/hcuts;

dels = {xdel ydel sdel hdel};
mins = {xmin ymin smin hmin};
maxs = {xmax ymax smax hmax};
%% Ps,Pyx definition
% Channel property Pxy and source Ps are defined
sig = 1;
iPs = [smin + sdel:sdel:smax];
Ps = cdf('norm',iPs,0,1) - [0 cdf('norm',[iPs(2:length(iPs))-sdel],0,1)];

iPx = [xmin + xdel:xdel:xmax]';
iPy = [ymin + ydel:ydel:ymax];
Pxy1 = cell2mat(arrayfun(@(x) (cdf('norm',iPy,x,1) - [0 cdf('norm',iPy(2:length(iPy))-ydel,x,1)]),iPx(1:ceil(length(iPx)/2)),'UniformOutput',false));
Pxy2 = cell2mat(arrayfun(@(x) (abs(1 - cdf('norm',iPy,x,1) - [1 - cdf('norm',[iPy(1:length(iPy) - 1)+ydel],x,1) 0])),iPx(ceil(length(iPx)/2) + 1:length(iPx)),'UniformOutput',false));
Pxy = cat(1,Pxy1,Pxy2);

%% Call to convexopt.m
% to compute Aeq,beq,A,b
[intmats,intmatC1,intmatC2,intmatC3] = convexopt(xdel,ydel,sdel,hdel,xi,yi,si,hi,Pxy);

%% definition of dels, ni
dels = [xdel ydel sdel hdel];
ni = [xi yi si hi];

%% Call to fmincon.m
Aeq = cat(1,intmatC1,intmatC2,intmatC3);
beq = cat(1,transpose(Ps),[1],zeros(xi*yi,1));
options = optimset('Algorithm','interior-point','Display','final-detailed','MaxFunEvals',1000000);
[lb x0] = deal(zeros(xi*yi*si*hi,1));
[x fval] = fmincon(@(Q) costfun(Q,dels,mins,maxs),x0,[],[],Aeq,beq,lb,[],@(Q) nonlcon(Ps,Pxy,Q,dels,ni,intmats) ,options);

% compact display
x = x';

% alert user
figure();
%% NONLCON
function [c,ceq] = nonlcon(Ps,Pxy,Q,dels,ni,intmats)
%% PARAMETERS
dcell = num2cell(dels);
[xdel,ydel,sdel,hdel] = dcell{:};

xi = ni(1);
yi = ni(2);
si = ni(3);
hi = ni(4);

li = xi*yi*si*hi;

[intmat1 intmat2 intmat3 intmat4 intmat5 intmat6] = intmats{:};

%% CONSTRAINT4 : MUTUAL INFORMATION CONSTRAINT
% I(aPs) - I(bP(y/x) ) <= 0 step1: a,b ==> integral form of Q ==> compute intmats
Qxyhs = rearrange(Q,xi,yi,si,hi);

xyQ = (xdel*ydel)*intmat2*(intmat1*Q);
xyQm = reshape(xyQ,si,hi);
syxQ = sum(xyQm,1);
% xyQ = (li/xi.yi , 1) syxQ = (1,hi) Ps = (1,si)
% xyQm = (si,hi)

repPs = repmat(transpose(Ps),1,hi);
Ahs = xyQm ./ repPs;
% Ahs = (si,hi) = a(h/s)

numPs = xyQm;
denomPs = transpose(Ps) * syxQ;
% (si,hi) = (si,1) * (1,hi)

IaPs = sdel*hdel*sum(sum(numPs.*log2(numPs./denomPs)));

% IbP
shQ = intmat5*(intmat4*transpose(Qxyhs));
shQm = reshape(shQ,yi,xi);
% shQm = (yi,xi)

hsQ = intmat5*(intmat4*transpose(Qxyhs));
hsQm = reshape(shQ,yi,xi);
yhsQ = (ydel*sdel*hdel)*sum(hsQm,1);
% yhsQ = (1,xi)

bx = transpose(yhsQ);
% bx = (xi,1)

repBx = repmat(bx,1,yi);
productP = Pxy .* repBx;

xProductP = xdel*sum(productP,1);
% (1,yi)
denomP = bx * xProductP;
% (xi,yi) = (xi,1) * (1,yi)

IbP = xdel*ydel*sum(sum(productP.*log2(productP./denomP)));
    
ineq1 = IaPs - IbP;

%% Ceq,C

ceq = [];
c = ineq1;

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

end % main

%% CONVEXOPT for linear constraints
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
toc;
end % witsenhausen
