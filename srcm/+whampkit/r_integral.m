function rc = r_integral(y, al)
%R_INTEGRAL Numerical integration of R_FUNCTION.
%   Two 16-point Gauss-Legendre quadratures -- over the finite interval UL
%   and over the exponential tail XO -- are combined and then corrected
%   analytically.  R is odd in Y up to the sign SIG, so a negative argument
%   is reflected and the sign restored at the end.
[abscissa, weight] = whampkit.gauss_legendre_16();
if real(y) < 0
    y = -y;
    sig = -1;
else
    sig = 1;
end
ya = imag(y);
yr = real(y);
ul = pi-real(2.8*y/(36+y));
coterm = cot(pi*y);
dv = pi*(1+coterm^2);          % d(cot(pi*y))/dy contribution
c = yr/al;
xo = log(c+sqrt(1+c^2));       % asinh(yr/al)
rc = complex(zeros(1,4));
for i = 1:16
    x = ul/2*(1+abscissa(i));
    z = sin(x);
    c = cos(x);
    g = yr/al*x/z;
    t = sqrt(1+g^2);
    b = log(g+t);
    g = (1/x-c/z)*g/t;
    t = al*(t*c-1);
    z = real(exp(x*ya));
    c = 0.5*(z+1/z);
    s = 0.5i*(z-1/z);
    f = coterm+g;
    h = 1-g*coterm;
    exf = exp(t-y*b);
    o = b*c+x*s;
    pv = x*c-b*s;
    xy = x*yr;
    r = (f*c+h*s)*exf;
    ry = (f*o-h*pv+dv*(c-g*s))*exf;
    rp = ((f*t-h*xy)*c+(h*t+f*xy)*s)*exf;
    rpy = (f*(t*o-xy*pv)-h*(t*pv+xy*o) ...
        +((t+xy*g)*c-(g*t-xy)*s)*dv)*exf;
    x = xo/2*(1+abscissa(i));
    z = exp(x);
    c = (z+1/z)/2-1;
    pv = exp(al*c-y*x);
    rc(1) = rc(1)+weight(i)*(ul*r+xo*pv);
    rc(2) = rc(2)-weight(i)*(ul*ry+xo*x*pv);
    rc(3) = rc(3)+weight(i)*(ul*rp+xo*al*c*pv);
    rc(4) = rc(4)-weight(i)*(ul*rpy+xo*al*x*c*pv);
end
% Analytic endpoint corrections, with the SIG reflection undone.
o = y/al;
pv = y^2/2;
rc(1) = o*(y*rc(1)/2-1)*sig;
rc(2) = 2*rc(1)+o*(pv*rc(2)+1)*sig;
rc(3) = y*o*rc(3)/2*sig;
rc(4) = 2*rc(3)+o*pv*rc(4)*sig;
end

