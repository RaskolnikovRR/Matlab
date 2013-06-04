function [x,fval] = optimiser()
%% function definition and intervals
[xmin,ymin,hmin,smin] = deal(0);
[xmax,smax,hmax,ymax] = deal(3);
[xcuts,ycuts,scuts,hcuts] = deal(6);

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

% xA = cat(2,xA,x);
% fvalA = cat(1,fvalA,fval');
% nonlcon(Ps,Pxy,x,dels,ni,intmats)
% 
% display('Q value at minimum obj fun (each column is for a different p):');
% display(xA);
% display('Value of obj fun at different p:');
% display(fvalA);