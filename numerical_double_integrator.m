function [intResult] = numerical_double_integrator(f,xmin,xmax,ymin,ymax,xcuts,ycuts)

syms i x1 x2 y1 y2 intResult;
% INSTRUCTION: input function f to be a character expression of x,y

% converting string f to sym expression f
f=sym(f);

% calculation of x and y update intervals
xdel = (xmax - xmin)/xcuts;
ydel = (ymax - ymin)/ycuts;

% start
x1 = xmin;
y1 = ymin;
x2 = Xupdate(x2,xdel);
y2 = Yupdate(y2,ydel);

intResult = 0;
i=1;
j=1;
while x2 <= xmax
    y1 = ymin;
    y2 = y1 + ydel;
    while y2 <= ymax
        intResult= intResult + addPiece(x2,y2,i,j);
        innerIterate();
    end
    outerIterate();
        % the current status is stored in a file rather than on terminal
        log=fopen('integrationlog','a');
        fprintf(log,'%s %s %s\n',char(vpa(x2)),char(vpa(y2)),char(vpa(intResult)));
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

    function piece = addPiece(x2,y2)
        piece = vpa(xdel*ydel*subs(f,{'x','y'},{x2 - xdel/2,y2 - ydel/2}),16);
    end

end

        
        