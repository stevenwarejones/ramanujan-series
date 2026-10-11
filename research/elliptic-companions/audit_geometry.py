#!/usr/bin/env python3
"""Exact SymPy geometry checks."""
import json
import sympy as S
x,y=S.symbols('x y'); i=S.I; r=1+S.sqrt(2); alpha=S.expand(r*r); beta=3-2*S.sqrt(2)
p=S.expand((1+x*x)*(alpha+beta*x*x))
X=S.sqrt(2)*(x+i)/(x-i)-1
Y=S.sqrt(2)*i*y/(x-i)**2
assert S.simplify((Y*Y-X**3+X).subs(y*y,p))==0
assert S.simplify(S.diff(X,x)/Y+2/y)==0
assert S.simplify(-i*Y/(X+1)-y/(1+x*x))==0
T=(alpha+beta*x*x)/(1+x*x)
assert S.simplify(T*T.subs(x,alpha/x)-1)==0
points=[(-r,-i*S.sqrt(2)*r), (i,-1+i)]
for point, expected in zip(points,[(-1,0),(0,0)]):
    xx,yy=point
    assert S.simplify(yy*yy-xx**3+xx)==0
    slope=(3*xx*xx-1)/(2*yy)
    dx=S.simplify(slope*slope-2*xx)
    dy=S.simplify(slope*(xx-dx)-yy)
    assert S.simplify(dx-expected[0])==0 and dy==expected[1]
for xx,yy,point in [(0,r,points[0]),(r,alpha+1,points[1])]:
    assert S.simplify(X.subs(x,xx)-point[0])==0
    assert S.simplify(Y.subs({x:xx,y:yy})-point[1])==0

import sympy as s
X,Y=s.symbols('X Y');i=s.I
G=-i*Y/(X+1)
rotated=G.subs({X:-X,Y:i*Y},simultaneous=True)
assert s.cancel((G*rotated-i*X).subs(Y**2,X**3-X))==0
# Rational coefficient of dX in dlog(g), using the curve equation.
eta=(3*X**2-1)/(2*(X**3-X))-1/(X+1)
xi=1/X
assert s.cancel(-eta.subs(X,-X)-(xi-eta))==0
assert s.cancel(-xi.subs(X,-X)-xi)==0
# Translation by (1,0) inverts g; translation by (-1,0) negates g.
for xx,yy,target in [((X+1)/(X-1),-2*Y/(X-1)**2,1/G),((1-X)/(X+1),-2*Y/(X+1)**2,-G)]:
    val=G.subs({X:xx,Y:yy},simultaneous=True)
    assert s.cancel(s.together(val-target).subs(Y**2,X**3-X))==0
# Differential independence follows already from unequal pole residues.
residues=s.Matrix([[1,0],[1,2],[-1,0],[-1,-2]])
assert residues.rank()==2
print(json.dumps({'curve':'Y^2=X^3-X','geometry_checks':'birational map, differential, logarithm, involution, endpoints and both exact order-four doublings passed','rotation':'(X,Y)->(-X,iY)',
'product_identity':'g(rho(P))*g(P)=iX(P)', 'differential_action':['rho^* eta=xi-eta','rho^* xi=xi'],
'letters':['eta=dlog(g)','xi=dlog(X)'],'residue_matrix_points':['(1,0)','(0,0)','(-1,0)','O'],
'residue_matrix':[[int(v) for v in row] for row in residues.tolist()], 'residue_rank':2,
'half_period_actions':['translation_(1,0): g->1/g','translation_(-1,0): g->-g'],
'scope':'Exact symbolic audit; rotation enlarges the one-letter logarithmic alphabet and changes the path. No special-value independence theorem.'},indent=2))
