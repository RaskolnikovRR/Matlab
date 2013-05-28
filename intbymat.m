function [iRes] = intbymat(f,xmin,xmax,ymin,ymax,xcuts,ycuts)

syms xdel ydel i j x y;

f = sym(f);

ydel = (ymax -ymin)/ycuts;
xdel = (xmax - xmin)/xcuts;

% compute F matrix

x1 = xmin;
x2 = xmin + xdel;
i = 1;
j = 1;

while x2 < xmax
    y1 = ymin;
    y2 = ymin + ydel;
    while y2 < ymax
        F(i,j) = subs(f,{'x','y'},{(x2+x1)/2,(y2+y1)/2});
        j = j + 1;
        y1 = y2;
        y2 = y1 + ydel;
    end
    i = i + 1;
    x1 = x2;
    x2 = x1 + xdel;
end

% compute integration result from integration matrix
iRes = sum(sum(F));

         