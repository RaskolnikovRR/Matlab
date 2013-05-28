function iRes = int4vars(f,xmin,xmax,ymin,ymax,smin,smax,hmin,hmax)

syms x y h s;
f = sym(f);

int1 = int(f,x,xmin,xmax);
int2 = int(int1,y,ymin,ymax);
int3 = int(int2,h,hmin,hmax);
iRes = int(int3,s,smin,smax);

end
