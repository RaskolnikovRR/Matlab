function [out] = signum(x)
syms i;
i = 1;
while i <= length(x)
	if( double(x(i)) > 0 )
		out(i) = 1;
	elseif( double(x(i)) < 0 )
		out(i) = -1;
	else
		out(i) = 0;
    end
    i=i+1;
end