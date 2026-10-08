"""Exact maps and norm substitutions, with primitive integer coefficients."""
from flint import fmpz_poly
X=fmpz_poly([0,1])
MAPS={2:(X**3,fmpz_poly([0,-65536,52992,-3456]),fmpz_poly([16777216,28311552,15925248,2985984])),
      3:(X**4,fmpz_poly([0,-1259712,1469664,-317952]),fmpz_poly([136048896,725594112,1289945088,764411904]))}
def transform(H,N):
 h=H.degree()
 if N==4:
  num=64*(4-X)**3;den=X*X;q=fmpz_poly([int(H[h])]);power=fmpz_poly([1])
  for i in range(h-1,-1,-1):power*=den;q=q*num+int(H[i])*power
  return q
 a,b,c=MAPS[N];u=fmpz_poly([]);v=fmpz_poly([1]);power=fmpz_poly([1])
 for i in range(h-1,-1,-1):
  power*=a;u,v=-b*u+a*v,-c*u+int(H[i])*power
 q=c*u*u-b*u*v+a*v*v;div=a**(h+1)
 out,rem=divmod(q,div);assert not rem
 return out
def bpoly(p,d):
 m=p.degree();q=fmpz_poly([int(p[m])]);inn=fmpz_poly([d,0,-1])
 for i in range(m-1,-1,-1):q=q*inn+int(p[i])*d**(m-i)
 return q//q.content()
def canon(p):
 p=p//p.content()
 return p if p[p.degree()]>0 else -p
