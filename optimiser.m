function [x,fval,elapse] = optimiser(gaps)
%% function definition and intervals
tic;
[xmin,ymin,hmin,smin] = deal(-6);
[xmax,smax,hmax,ymax] = deal(6);
[xcuts,ycuts,scuts,hcuts] = deal(gaps);

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
cvar = 1;
sig = 2;
iPs = [smin + sdel:sdel:smax];
Ps = cdf('norm',iPs,0,sig) - [0 cdf('norm',[iPs(2:length(iPs))-sdel],0,sig)];

iPx = [xmin + xdel:xdel:xmax]';
iPy = [floor(ymin + ydel - 4*cvar):ydel:ceil(ymax + 4*cvar)];
Pxy = cell2mat(arrayfun(@(x) (cdf('norm',iPy,x,1) - [0 cdf('norm',iPy(2:length(iPy))-ydel,x,1)]),iPx,'UniformOutput',false));
Pxy(Pxy == 0) = realmin;
ysig = size(Pxy,2);
yi = ysig;

%% Call to convexopt.m
% to compute Aeq,beq,A,b
[intmats,intmatC1,intmatC2,intmatC3] = convexopt(ysig,xdel,ydel,sdel,hdel,xi,yi,si,hi,Pxy);

%% definition of dels, ni
dels = [xdel ydel sdel hdel];
ni = [xi yi si hi];

%% Call to fmincon.m
Aeq = cat(1,intmatC1,intmatC2,intmatC3);
beq = cat(1,transpose(Ps),[1],zeros(xi*ysig,1));
options = optimset('Algorithm','interior-point','Display','final-detailed','MaxFunEvals',1000000,'MaxIter',4000);
[lb] = deal(zeros(xi*yi*si*hi,1));
[x0] = 0.016*ones(xi*yi*si*hi,1);
[x fval] = fmincon(@(Q) costfun(Q,dels,mins,maxs),x0,[],[],Aeq,beq,lb,[],@(Q) nonlcon(ysig,Ps,Pxy,Q,dels,ni,intmats) ,options);
elapse = toc; % elapsed time
figure(); % user alert
title('SUCCESS');
end % main