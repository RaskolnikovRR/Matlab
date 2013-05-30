function cost = costfun(Q)
%% COST FUNCTIN definition
function k = kk(s,h)
if s == h
    k = 0;
else 
    k = 1;
end
end
    
%% Loop to compute kappa
hi = 1;
i = 1;

while( hi <= 2)
	si = 1;
	while( si <= 2)
		yi = 1;
		while( yi <= 2)
			xi = 1;
			while( xi <= 2)
				intg(i) = kk(si,hi)*Q(i);
				i = i + 1;
				xi = xi + 1;
			end
			yi = yi + 1;
		end
		si = si + 1;
	end
	hi = hi + 1;
end

%% INTEGRATION OVER s,x,y,h
cost = sum(intg);

end % main
