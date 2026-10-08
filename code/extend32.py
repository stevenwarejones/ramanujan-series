"""Search for counterexamples to the N=73117 degree-32 candidate.

For h>D and deg(x)<=D, j cannot lie in Q(x). Since j satisfies a quadratic
over Q(x), the remainder H(j)=(u(x)j+v(x))/a(x)^h forces u(x)=v(x)=0.
Thus gcd(u,v), rather than the full norm resultant, contains every eligible x.
No conjecture about factor degrees or isogeny splitting is needed for this cut.
"""
from run_paths import output_path,generated_input
import json,time,math
from pathlib import Path
from flint import fmpz_poly,arb,ctx
from polynomials import MAPS,transform,bpoly,canon
ROOT=Path(__file__).resolve().parents[1]

def low_degree_cover(H,N,D):
    h=H.degree()
    if N==1:
        return canon(fmpz_poly([int(H[h-i])*1728**(h-i) for i in range(h+1)]))
    if h<=D or N==4:return canon(transform(H,N))
    a,b,c=MAPS[N];u=fmpz_poly([]);v=fmpz_poly([1]);power=fmpz_poly([1])
    for i in range(h-1,-1,-1):
        power*=a;u,v=-b*u+a*v,-c*u+int(H[i])*power
    return canon(u.gcd(v))

def remove_zero(p):
    co=list(p.coeffs())
    while co and co[0]==0:co.pop(0)
    return fmpz_poly(co)

def rouche(p,den):
    m=p.degree()
    return abs(int(p[0]))*den**m > sum(abs(int(p[i]))*den**(m-i) for i in range(1,m+1))

def main():
    ctx.dps=100;start=time.perf_counter();D=32;den=10**366
    cat=json.loads(generated_input('catalogue64.json').read_text())
    can=json.loads((ROOT/'results'/'candidate32.json').read_text());winner=fmpz_poly(can['p'])
    roots=winner.complex_roots();wz=min((z for z,e in roots),key=lambda z:float(abs(z).log()))
    assert abs(wz)<arb(1)/den
    assert all(abs(wz)<abs(z) for z,e in roots if not z.overlaps(wz))
    excluded=0;factor_sets=0;empty=0;potential=[];ledger=[]
    for index,o in enumerate(reversed(cat['orders'])):
        d=o['d'];h=o['h'];H=fmpz_poly.hilbert_class_poly(-d);assert H.degree()==h
        row={'d':d,'h':h,'families':{}}
        for N in (1,2,3,4):
            if N in (1,4) and h>D:continue
            q=remove_zero(low_degree_cover(H,N,D))
            if q.degree()<=0:
                empty+=1;row['families'][str(N)]='empty-cover';continue
            if rouche(q,den):
                excluded+=1;row['families'][str(N)]='integer-Rouche';continue
            factor_sets+=1;content,fs=q.factor();check=fmpz_poly([content])
            for p,e in fs:check*=p**e
            assert check==q
            dispositions=[]
            for p,e in fs:
                p=canon(p);m=p.degree()
                if m>D:dispositions.append([m,'degree']);continue
                if m==1 and p(1)==0:dispositions.append([m,'boundary']);continue
                if rouche(p,den):dispositions.append([m,'Rouche']);continue
                bp=bpoly(p,d);bc,bfs=bp.factor();bcheck=fmpz_poly([bc])
                for f,ee in bfs:bcheck*=f**ee
                assert bcheck==bp
                degs=sorted((f.degree(),ee) for f,ee in bfs)
                assert degs in [[(2*m,1)],[(m,1),(m,1)]],(N,d,degs)
                bd=degs[0][0]
                if bd>D:dispositions.append([m,'B-degree',bd]);continue
                zs=p.complex_roots()
                if list(map(int,p.coeffs()))==can['p']:
                    dispositions.append([m,'equals-candidate',bd]);continue
                for z,ee in zs:
                    if abs(wz)<abs(z):continue
                    potential.append(dict(level=N,d=d,h=h,degree=m,full_x_B_degree=bd,x=str(z),rate=str(-abs(z).log()/arb(10).log()),p=list(map(int,p.coeffs()))))
                dispositions.append([m,'balls',bd])
            row['families'][str(N)]=dispositions
        ledger.append(row)
        if index%500==0:
            print('EXTEND32',index,'/',len(cat['orders']),'d',d,'seconds',round(time.perf_counter()-start,1),'factor sets',factor_sets,'challengers',len(potential),flush=True)
    result=dict(status='PASS' if not potential else 'CHALLENGERS',budget=32,orders=len(cat['orders']),
                empty_covers=empty,Rouche_exclusions=excluded,factor_sets=factor_sets,
                candidate_N=73117,potential=potential,ledger=ledger,elapsed=time.perf_counter()-start)
    output_path('extension32.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='ledger'}),flush=True)

if __name__=='__main__':main()
