function[iResult] = numInt(func)
syms a b limits del resultArray i;
limits = 10;
del = 1;
i = 1;
func = str2func(func);
while(1) % manual interrupt when 'satisfiable' results are achieved
	iResult = 0;
	a = limits*(-1)
	input('Continue iterations?');
	while( a < limits)
		b = a + del;
		iResult = iResult + func( (a+b)/2) * del;
		a = b;
	end
	resultArray(i) = iResult
	dlmwrite('Numerical Integration iterations','-append');
	i = i + 1;
	limits = limits*10;
end	
	