function rc = r_asymptotic(y, al)
%R_ASYMPTOTIC Asymptotic expansion of R_FUNCTION for large AL, moderate |Y|.
%   Sums the P and Q series until the terms start growing again.
y2 = y^2;
coterm = cot(pi*y);
pn = -y/al;
pyn = pn;
a = 1/(al*sqrt(2*pi*al));
qn = pi*y2*coterm*a;
qyn = qn*(2-y*pi*coterm)-y*pi^2*y2*a;
p = pn+qn;
py = pyn+qyn;
pp = -pn-1.5*qn;
ppy = -pyn-1.5*qyn;
ay = abs(y)+2;
t = 0;
c = 7.2e35;
for n = 1:100
    m = n-1;
    pyn = (pyn*(m^2-y2)-2*y2*pn)/((2*m+1)*al);
    pn = pn*(m^2-y2)/((2*m+1)*al);
    qyn = (qyn*((m+0.5)^2-y2)-2*y2*qn)/(2*n*al);
    qn = qn*((m+0.5)^2-y2)/(2*n*al);
    if m >= ay
        c = n*(abs(pn)+abs(qn));
        if c <= 1e-7*abs(pp) || c >= t
            break
        end
    end
    p = p+pn+qn;
    py = py+pyn+qyn;
    pp = pp-(n+1)*pn-(n+1.5)*qn;
    ppy = ppy-(n+1)*pyn-(n+1.5)*qyn;
    t = c;
end
rc = [p+pn+qn, py+pyn+qyn, pp+p, ppy+py];
end

