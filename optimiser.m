function [points,results] = optimiser()

options = optimset('Algorithm','interior-point');

p = 0.1;
i = 1;
ilims = 10;
while (i < ilims && p < 1)
    [points(i,:) results(i)] = fmincon(@(Q) costfun(Q),zeros(1,16),[],[],[],[],zeros(1,16),[],@(Q) nonlcon(p,Q),options);
    p = p + 1/ilims;
    i = i + 1;
end
