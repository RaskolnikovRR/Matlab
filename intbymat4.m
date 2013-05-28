function [xi,yi,si,hi,iRes,F] = intbymat4(f,xmin,xmax,ymin,ymax,smin,smax,hmin,hmax,xcuts,ycuts,scuts,hcuts)

syms xdel ydel sdel hdel i;
syms h s x y;
syms hi si xi yi;

ydel = (ymax -ymin)/ycuts;
xdel = (xmax - xmin)/xcuts;
sdel = (smax -smin)/scuts;
hdel = (hmax - hmin)/hcuts;

i = 1;
[hi,si,xi,yi] = deal(0);

% initialization
h1 = hmin; h2 = h1 + hdel;
s1 = smin; s2 = s1 + sdel;
x1 = xmin; x2 = x1 + xdel;
y1 = ymin; y2 = y1 + ydel;

while h2 <= hmax
    h = avg(h1,h2);
    s1 = smin; s2 = s1 + sdel; si = 0;
    while s2 <= smax
        s = avg(s1,s2);
        y1 = ymin; y2 = y1 + ydel; yi = 0;
        while y2 <= ymax
            y = avg(y1,y2);
            x1 = xmin; x2 = x1 + xdel; xi = 0;
            while x2 <= xmax
                x = avg(x1,x2);
                F(i) = f(x,y,s,h);
                i = i + 1;
                xi = xi + 1;
                [x1 x2] = varupdate(xdel,x2);
            end
            % the current status is stored in a file rather than on terminal
            log=fopen('intbymat4log.txt','a');
            fprintf(log,'%d %d %d %d %d\n',x2,y2,h2,s2,i);
            fclose(log);
            [y1 y2] = varupdate(ydel,y2);
            yi = yi +1;
        end
        log=fopen('intbymat4log.txt','a');
        fprintf(log,'%d %d %d %d %d\n\n',x2,y2,h2,s2,i);
        fclose(log);
        [s1 s2] = varupdate(sdel,s2);
        si = si + 1;
    end
    log=fopen('intbymat4log.txt','a');
    fprintf(log,'%d %d %d %d %d\n\n',x2,y2,h2,s2,i);
    fclose(log);
    [h1 h2] = varupdate(hdel,h2);
    hi = hi + 1;
end

length(F)
iRes = (xdel*ydel*hdel*sdel)*sum(F);

    function [l u] = varupdate(del,u)
        l = u;
        u = l + del;
    end

    function m = avg(a,b)
        m = (a+b)/2;
    end

end