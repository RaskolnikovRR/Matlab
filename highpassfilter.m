syms t X fre;
t = linspace(0,10,0.1);
X = sin ( 2*pi*t);
fre = fft(X);
fre