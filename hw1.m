clear;clc;close all;

syms x;
pu = 4*10^4*(x-1)^2+5.4*10^4;
pl = 2*10^4*(x-1)^2+1.73*10^5;
tu = 288*x^(-0.2);
tl = 731*x^(-0.2);
th = 0;
LE = 0;
TE = 1;
a = 10;
yu = 0;
yl = 0;

N = Nprime(LE,TE,pu,pl,tu,tl,th);
N = sym2poly(N);

A = Aprime(LE,TE,pu,pl,tu,tl,th);
A = sym2poly(A);

L = N*cosd(a) - A*sind(a);
D = N*sind(a) + A*cosd(a);

Mle = Mleprime(LE,TE,pu,pl,tu,tl,th,yu,yl);
Mle = sym2poly(Mle);

xcp = -Mle/N;

Ml4 = Mle + (TE-LE)/4 * L;

function N = Nprime(LE,TE,pu,pl,tu,tl,th)

    expr1 = (pu*cosd(th) + tu*sind(th));
    expr2 = (pl*cosd(th) - tl*sind(th));
    
    N = -1 * (int(expr1,LE,TE)) + int(expr2,LE,TE);

end

function A = Aprime(LE,TE,pu,pl,tu,tl,th)

    expr1 = (-pu*sind(th) + tu*cosd(th));
    expr2 = (pl*sind(th) + tl*cosd(th));
    
    A = int(expr1,LE,TE) + int(expr2,LE,TE);

end
function M = Mleprime(LE,TE,pu,pl,tu,tl,th,yu,yl)

    syms x;

    expr1 = (pu*cosd(th) + tu*sind(th))*x - (pu*sind(th) - tu*cosd(th))*yu;
    expr2 = (-pl*cosd(th) + tl*sind(th))*x + (pl*sind(th) + tl*cosd(th))*yl;
    
    M = int(expr1,LE,TE) + int(expr2,LE,TE);

end