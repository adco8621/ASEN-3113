function M = Mleprime(LE,TE,pu,pl,tu,tl,th,yu,yl)

    syms x;

    expr1 = (pu*cosd(th) + tu*sind(th))*x - (pu*sind(th) - tu*cosd(th))*yu;
    expr2 = (-pl*cosd(th) + tl*sind(th))*x + (pl*sind(th) + tl*cosd(th))*yl;
    
    M = int(expr1,LE,TE) + int(expr2,LE,TE);

end