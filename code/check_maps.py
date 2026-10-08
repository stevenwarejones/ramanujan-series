"""Prove the polynomial j/x maps directly by substitution in rational t maps."""
from flint import fmpz_poly
from polynomials import MAPS

def homogenize(p,n,d,m):
    return sum((int(p[i])*n**i*d**(m-i) for i in range(p.degree()+1)),fmpz_poly([]))

T=fmpz_poly([0,1])
for N,c in ((2,64),(3,27)):
    xn=4*c*T;xd=(T+c)**2
    jn=(T+256)**3 if N==2 else (T+27)*(T+243)**3
    jd=T**N
    a,b,cc=MAPS[N];m=max(v.degree() for v in (a,b,cc))
    assert homogenize(a,xn,xd,m)*jn**2+homogenize(b,xn,xd,m)*jn*jd+homogenize(cc,xn,xd,m)*jd**2==0
    # x(c^2/t)=x(t), clearing the respective denominators.
    assert (4*c*c*c*T)*(T+c)**2 == (4*c*T)*(c*c+c*T)**2
x=4*T*(1-T)
assert 64*(4-x)**3*(T*T*(1-T)**2)==256*(1-T+T*T)**3*x*x
print('Exact modular rational maps PASS')
