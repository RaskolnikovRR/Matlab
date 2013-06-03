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
