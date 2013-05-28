function [iRes,ress] = quadnumint(xmin,xmax,ymin,ymax,xcuts,ycuts,s,k)

f = @(x,v) (exp(-x^2/(2*s^2))*exp(-v^2/2)/(2*pi*s)*(k^2*(s*signum(x) - x)^2 + ( s*signum(x) - s*tanh( s*(s*signum(x) + v)))^2 ) );
syms i x1 x2 y1 y2 iRes ress;
xdel = (xmax - xmin)/xcuts
ydel = (ymax - ymin)/ycuts
input('Press enter to continue:\n');
% f = str2func(f);
x1 = xmin;
y1 = ymin;
x2 = Xupdate(x2,xdel);
y2 = Yupdate(y2,ydel);

iRes = 0;
i=1;
j=1;
while x2 <= xmax
    y1 = ymin;
    y2 = y1 + ydel;
    while y2 <= ymax
        iRes= iRes + addPiece(x2,y2,i,j);
        innerIterate();
    end
    outerIterate();
        log=fopen('integrationlog17M13','a');
        fprintf(log,'%s %s %s\n',char(vpa(x2)),char(vpa(y2)),char(vpa(iRes)));
        fclose(log);
end

    function x2 = Xupdate(x2,xdel)
        x2 = x1 + xdel;
    end

    function y2 = Yupdate(y2,ydel)
        y2 = y1 + ydel;
    end

    function innerIterate()
        y1 = y2;
        y2 = y1 + ydel;
        i=i+1;
    end

    function outerIterate()
        x1 = x2;
        x2 = x1 + xdel;
        j=j+1;
    end

    function piece = addPiece(x2,y2,i,j)
        piece = vpa(xdel*ydel*f(x2 -xdel/2 ,y2-ydel/2 ),16);
    end

end

        
        