function rc = r_taylor(y, al)
%R_TAYLOR Taylor expansion of R_FUNCTION for small AL or large |Y|.
%   Summation stops once the term is eight orders below the leading one.
y2 = y^2;
pn = y/(y2-1);
pyn = -y*(y2+1)/(y2-1)^2;
rc = [pn pyn pn pyn];
for i = 2:100
    coterm = (2*i-1)/(y2-i^2)*al;
    pyn = coterm*(pyn-2*y2/(y2-i^2)*pn);
    pn = coterm*pn;
    rc = rc + [pn pyn i*pn i*pyn];
    if abs(pn)*1e8 < abs(rc(1))
        break
    end
end
end

