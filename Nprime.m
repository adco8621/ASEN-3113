function N = Nprime(LE,TE,pu,pl,tu,tl,th)

    expr1 = (pu*cosd(th) + tu*sind(th));
    expr2 = (pl*cosd(th) - tl*sind(th));
    
    N = -1 * (int(expr1,LE,TE)) + int(expr2,LE,TE);

end