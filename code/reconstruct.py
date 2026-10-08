"""Recompute every eligible factor, including candidates the discovery run discarded.

This is a second run with an independent conductor enumeration, not a second CAS.
FLINT factorization and certified Hilbert class polynomials are trusted dependencies.
"""
import json, time, hashlib
from pathlib import Path
from flint import fmpz_poly
from catalogue import build
from polynomials import transform, bpoly, canon
ROOT=Path(__file__).resolve().parents[1]

def key(r):return r['level'],r['d'],tuple(r['p'])
def checked_factor(p):
    content,fs=p.factor();q=fmpz_poly([content])
    for f,e in fs:q*=f**e
    assert q==p
    return fs

def main():
    start=time.perf_counter()
    cat=build(32)
    oldorders=json.loads((ROOT/'data'/'orders32.json').read_text())['orders']
    assert sorted(cat['orders'],key=lambda r:r['d'])==sorted(oldorders,key=lambda r:r['d'])
    old=json.loads((ROOT/'data'/'search_results.json').read_text())
    expected={key(r):r for r in old['records']};found={};ledger=[]
    old1=json.loads((ROOT/'data'/'level1_frontier.json').read_text())
    expected1={r['d']:r for r in old1['all_orders']};seen1=set()
    for index,o in enumerate(cat['orders']):
        d=o['d'];h=o['h'];H=fmpz_poly.hilbert_class_poly(-d)
        assert H.degree()==h
        hr=dict(d=d,h=h,H_sha256=hashlib.sha256(str(H).encode()).hexdigest(),families={})
        if h<=16 and d not in (3,4):
            p=canon(fmpz_poly([int(H[h-i])*1728**(h-i) for i in range(h+1)]))
            fs=checked_factor(bpoly(p,d));degs=sorted((f.degree(),e) for f,e in fs)
            assert degs in [[(2*h,1)],[(h,1),(h,1)]]
            assert expected1[d]['B_degree']==degs[0][0]
            seen1.add(d);hr['families']['1']={'parameter_degree':h,'B_degree':degs[0][0]}
        for N in (2,3,4):
            if N==4 and h>16:continue
            q=canon(transform(H,N));fs=checked_factor(q);entries=[]
            for p,e in fs:
                p=canon(p);m=p.degree();entry={'degree':m,'multiplicity':e}
                if not p[0] or (m==1 and p(1)==0):entry['disposition']='zero-or-boundary'
                elif m>16:entry['disposition']='parameter-degree-exceeds-16'
                else:
                    bfs=checked_factor(bpoly(p,d));degs=sorted((f.degree(),ee) for f,ee in bfs)
                    assert degs in [[(2*m,1)],[(m,1),(m,1)]]
                    bd=degs[0][0];entry['B_degree']=bd
                    if bd>16:entry['disposition']='joint-x-B-degree-exceeds-16'
                    else:
                        r=dict(level=N,d=d,h=h,parameter_degree=m,B_degree=bd,p=list(map(int,p.coeffs())))
                        k=key(r);assert k not in found;found[k]=r
                        assert k in expected,(N,d,m,bd,'missing from discovery')
                        assert expected[k]['B_degree']==bd
                        entry['disposition']='eligible'
                entries.append(entry)
            hr['families'][str(N)]=entries
        ledger.append(hr)
        if index%250==0:
            print('RECONSTRUCT',index,'/',len(cat['orders']),'eligible',len(found),'seconds',round(time.perf_counter()-start,1),flush=True)
    assert set(found)==set(expected)
    assert seen1==set(expected1)
    result=dict(status='PASS',orders=len(cat['orders']),fundamental_fields=cat['fields'],
                eligible_levels_2_3_4=len(found),checked_level1_orders=len(seen1),
                regenerated_every_candidate=True,checked_all_factor_products=True,
                elapsed=time.perf_counter()-start,ledger=ledger)
    (ROOT/'results'/'reconstruction.json').write_text(json.dumps(result,indent=2))
    print('PASS',json.dumps({k:v for k,v in result.items() if k!='ledger'}),flush=True)

if __name__=='__main__':main()
