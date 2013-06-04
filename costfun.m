function cost = costfun(Q,dels,mins,maxs)
%% cost function definition
kappa = @(x,y,s,h) ((s-h)^2+x^2);

%% Loop to compute kappa
dels = num2cell(dels);

[xdel ydel sdel hdel] = dels{:};
[xmin ymin smin hmin] = mins{:};
[xmax ymax smax hmax] = maxs{:};

hi = hmin + hdel;
i = 1;

while( hi <= hmax)
	si = smin + sdel; 
	while( si <= smax)
		yi = ymin + ydel;
		while( yi <= ymax)
			xi = xmin + xdel;
			while( xi <= xmax)
				intg(i) = kappa(xi,yi,si,hi)*Q(i);
				i = i + 1;
				xi = xi + xdel;
			end
			yi = yi + ydel;
		end
		si = si + sdel;
	end
	hi = hi + hdel;
end

%% INTEGRATION OVER s,x,y,h
cost = sum(intg);

end % main
