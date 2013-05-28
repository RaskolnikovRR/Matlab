function iRes = genmatrix(F,xi,yi,si,hi)

syms rowiter coliter li;

li = xi*yi*si*hi;

    function vec = intX(F)
        rowiter = 1;
        intmat1 = zeros(li/xi,li);
        
        while( rowiter <= li/xi )
            intmat1(rowiter,(xi*(rowiter-1)+1):xi*rowiter)=1;
            rowiter = rowiter + 1;
        end
        
        vec = intmat1*F;
    end

    function vec = intY(F)
        rowiter = 1;
        intmat2 = zeros(li/yi,li);
        
        while(rowiter <= li/yi)
            coliter = 1;
            while( coliter <= yi )
                intmat2(rowiter,coliter) = 1;
                coliter = nextY(coliter);
            end
        end
        vec = intmat2*F;
    end

    function vec = intS(F)
        rowiter = 1;
        intmat3 = zeros(li/si,li);
        
        while(rowiter <= li/si)
            coliter = 1;
            while( coliter <= si )
                intmat3(rowiter,coliter) = 1;
                coliter = nextS(coliter);
            end
        end
        
        vec = intmat3*F;
    end

    function vec = intH(F)
        intmat4 = zeros(li/hi,li);
        coliter = 1;
        
        while( rowiter <= li/hi)
            coliter = 1;
            while(coliter <= hi)
                intmat4(rowiter,coliter) = 1;
                coliter = nextH(coliter);
            end
        end
        
        vec = intmat4*F;
    end

    function i = nextH(i)
        i = i + hi;
    end

    function i = nextS(i)
        i = i + si;
    end

    function i = nextY(i)
        i = i + xi;
    end