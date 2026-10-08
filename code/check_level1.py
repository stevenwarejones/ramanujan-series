"""Exclude every level-1 conjugate against the global degree-16 frontier."""
import json
from pathlib import Path
from flint import fmpz_poly,arb,ctx
from polynomials import canon
from extend32 import rouche
ROOT=Path(__file__).resolve().parents[1];ctx.dps=120
one=json.loads((ROOT/'data'/'level1_frontier.json').read_text())
two=json.loads((ROOT/'data'/'search_results.json').read_text())['frontiers']['2']
wins={};checks=0;rootsets=0
for k in range(1,17):
    zs=fmpz_poly(two[k]['p']).complex_roots()
    w=min((z for z,e in zs),key=lambda z:float(abs(z).log()))
    assert all(abs(w)<abs(z) for z,e in zs if not w.overlaps(z))
    r=-abs(w).log()/arb(10).log();dec=int(float(r));assert dec<r<dec+1
    wins[k]=(w,10**dec)
for o in one['all_orders']:
    k=o['B_degree']
    if k>16:continue
    H=fmpz_poly.hilbert_class_poly(-o['d']);h=H.degree()
    p=canon(fmpz_poly([int(H[h-i])*1728**(h-i) for i in range(h+1)]))
    w,den=wins[k];assert abs(w)<arb(1)/den
    if not rouche(p,den):
        rootsets+=1
        assert all(abs(w)<abs(z) for z,e in p.complex_roots())
    checks+=1
print('Global frontier beats every level-1 conjugate:',checks,'orders;',rootsets,'root sets')
