function iRes = genmatrix(F,xi,yi,si,hi)

syms rowiter coliter li;

li = xi*yi*si*hi;

% initializations
rowiter = 1;
intmat1 = zeros(yi*si*hi,yi*hi*si*xi);

while( rowiter <= li/xi )
    intmat1(rowiter,(xi*(rowiter-1)+1):xi*rowiter)=1;
    rowiter = rowiter + 1;
end

rowiter = 1;
intmat2 = zeros(li/(xi*yi),li/xi);

while(rowiter <= li/(xi*yi))
    intmat2(rowiter,(yi*(rowiter-1)+1):yi*rowiter)=1;
    rowiter = rowiter + 1;
end

rowiter = 1;
intmat3 = zeros(li/(xi*yi*si),li/(xi*yi));

while(rowiter <= li/(xi*yi*si))
    intmat3(rowiter,(si*(rowiter-1)+1):si*rowiter)=1;
    rowiter = rowiter + 1;
end

iRes = (intmat3*intmat2*intmat1)*transpose(F);

end
