% turbine inlet temperature max = 620 celcius

T = [6,7,8,9,10,20,25,30,35];

for i = 1:numel(T)
    data(i,1)=T(i);
end
for i = 1:numel(T)
    fprintf('\niteration number: %d',i);
    [data(i,2),data(i,3)]= Regenerative(585,0.35,3,T(i),110);
    fprintf('---------------------------------------\n');
end

disp(data);