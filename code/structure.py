"""Exact checks of the ramification/genus-character mechanism and winner metadata."""
import json,math,time
from pathlib import Path
from flint import fmpz_poly,arb,acb,ctx
from catalogue import kronecker_prime,reduced_forms
from polynomials import transform,bpoly,canon
ROOT=Path(__file__).resolve().parents[1]

def main():
    start=time.perf_counter();ctx.dps=100
    old=json.loads((ROOT/'data'/'search_results.json').read_text())
    cat=json.loads((ROOT/'results'/'catalogue64.json').read_text())
    orders={o['d']:o for o in cat['orders']}
    winners=[]
    for k in range(1,17):
        r=old['frontiers']['2'][k];o=orders[r['d']]
        winners.append(dict(budget=k,N=r['d']//4,degree=r['B_degree'],h=r['h'],
                            fundamental_d=o['fundamental_d'],conductor=o['conductor'],rate=r['rate_mid']))
    Ns=sorted({r['N'] for r in winners}|{5,9,13,17,25,33,41,49,57,65,73,81,89,97})
    checks=[]
    for N in Ns:
        H=fmpz_poly.hilbert_class_poly(-4*N);h=H.degree()
        assert h==len(reduced_forms(4*N))
        q=canon(transform(H,2));fs=q.factor()[1]
        # Unique horizontal factor is repeated twice; select the actual tau branch.
        tau=acb(arb(-1)/2,arb(N).sqrt()/2)
        t=(tau.modular_eta()/(2*tau).modular_eta())**24;x=256*t/(t+64)**2
        hits=[]
        for p,e in fs:
            value=acb(0)
            for c in reversed(p.coeffs()):value=value*x+int(c)
            if value.contains(0):hits.append((canon(p),e))
        assert len(hits)==1
        p,e=hits[0];assert e==2 and 2*p.degree()==h
        for k in range(1,17):
            legacy=old['frontiers']['2'][k]
            if legacy['d']!=4*N:continue
            assert list(map(int,p.coeffs()))==legacy['p']
            zs=p.complex_roots();small=min((z for z,ee in zs),key=lambda z:float(abs(z).log()))
            assert small.overlaps(x)
            assert all(abs(small)<abs(z) for z,ee in zs if not z.overlaps(small))
        bfs=bpoly(p,N).factor()[1];bd=min(f.degree() for f,ee in bfs)
        predicted=h//2 if N%8==5 else h
        assert bd==predicted
        checks.append(dict(N=N,h=h,parameter_degree=p.degree(),full_x_B_degree=bd,
                           genus_sign=kronecker_prime(N,2),N_mod_8=N%8))
    result=dict(status='PASS',winner_rows=winners,character_checks=checks,
                exact_checks=len(checks),elapsed=time.perf_counter()-start,
                interpretation='Known CM/Fricke/genus theory explains both degree halving and radical absorption.')
    (ROOT/'results'/'structure.json').write_text(json.dumps(result,indent=2))
    print('STRUCTURE PASS',len(checks),'checks',time.perf_counter()-start)

if __name__=='__main__':main()
