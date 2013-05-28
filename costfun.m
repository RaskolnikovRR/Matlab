function cost = costfun(Q)

%~~~ for binary variables - kappa~~
k = @(s,x,h,y) ((s-h)^2+x^2);

[hi,si,xi] = deal(1);
i = 1;

while( hi <= 2)
	si = 1;
	while( si <= 2)
		yi = 1;
		while( yi <= 2)
			xi = 1;
			while( xi <= 2)
				intg(i) = k(si,xi,hi,yi)*Q(i);
				i = i + 1;
				xi = xi + 1;
			end
			yi = yi + 1;
		end
		si = si + 1;
	end
	hi = hi + 1;
end

% INTEGRATION OVER s,x,y,h
cost = sum(intg);
