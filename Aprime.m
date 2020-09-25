function A = Aprime(LE,TE,pu,pl,tu,tl,th)

    expr1 = (-pu*sind(th) + tu*cosd(th));
    expr2 = (pl*sind(th) + tl*cosd(th));
    
    A = int(expr1,LE,TE) + int(expr2,LE,TE);

end