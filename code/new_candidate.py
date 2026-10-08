"""Exact x/B certificate for the N=253*17^2 candidate.

A/B containment uses the separately stated CM descent lemma. The output is
not a proof of global degree-32 optimality and does not give an explicit A.
"""
from run_paths import output_path
import json,time,sys
from pathlib import Path
from flint import fmpz_poly,arb,acb,ctx
from polynomials import transform,bpoly,canon
ROOT=Path(__file__).resolve().parents[1]

def compose(p,q):
    r=fmpz_poly([])
    for c in reversed(p.coeffs()):r=r*q+c
    return r

def radical_certificate(p,N):
    bp=bpoly(p,N);content,fs=bp.factor();check=fmpz_poly([content])
    for f,e in fs:check*=f**e
    assert check==bp
    degs=sorted((f.degree(),e) for f,e in fs)
    if degs==[(p.degree(),1),(p.degree(),1)]:
        g=fs[0][0];E=fmpz_poly(list(g.coeffs())[::2]);O=fmpz_poly(list(g.coeffs())[1::2])
        r=fmpz_poly([N,-N]);num=-compose(E,r);den=compose(O,r)
        num%=p;den%=p
        assert den and (num*num-r*den*den)%p==0
        return dict(full_x_B_degree=p.degree(),B_factors=[list(map(int,f.coeffs())) for f,e in fs],
                    numerator=list(map(int,num.coeffs())),denominator=list(map(int,den.coeffs())),
                    identity='numerator^2 = N*(1-x)*denominator^2 modulo p; choose sign B>0')
    assert degs==[(2*p.degree(),1)]
    return dict(full_x_B_degree=2*p.degree(),B_factors=[list(map(int,fs[0][0].coeffs()))])

def main():
    start=time.perf_counter();N=73117;d=4*N;ctx.dps=460
    H=fmpz_poly.hilbert_class_poly(-d);assert H.degree()==64
    q=canon(transform(H,2));content,fs=q.factor();check=fmpz_poly([content])
    for f,e in fs:check*=f**e
    assert check==q
    tau=acb(arb(-1)/2,arb(N).sqrt()/2)
    t=(tau.modular_eta()/(2*tau).modular_eta())**24
    x=256*t/(t+64)**2
    assert x.imag.contains(0) and x.real<0
    matches=[]
    for p,e in fs:
        z=acb(0)
        for c in reversed(p.coeffs()):z=z*x+int(c)
        if z.contains(0):matches.append(canon(p))
    assert len(matches)==1
    p=matches[0];assert p.degree()==32
    cert=radical_certificate(p,N);assert cert['full_x_B_degree']==32
    # Distinguish the selected embedding using certified isolation.
    roots=p.complex_roots();hits=[z for z,e in roots if z.overlaps(x)]
    assert len(hits)==1 and hits[0].imag.contains(0)
    z=hits[0];assert all(abs(z)<abs(w) for w,e in roots if not w.overlaps(z))
    rate=-abs(z).log()/arb(10).log()
    print('NEW CANDIDATE',N,'h',H.degree(),'factor degrees',[(f.degree(),e) for f,e in fs],flush=True)
    print('rate',rate,flush=True)
    result=dict(N=N,d=d,h=64,conductor=17,fundamental_d=1012,parameter_degree=32,
                x=str(z),rate=str(rate),p=list(map(int,p.coeffs())),radical=cert,
                full_degree_status='32, conditional on the stated standard CM A/B descent lemma',
                optimality_status='not certified across all four families',
                explicit_A_status='not constructed in this certificate',elapsed=time.perf_counter()-start)
    output_path('candidate32.json').write_text(json.dumps(result,indent=2))
    print('saved',time.perf_counter()-start,flush=True)

if __name__=='__main__':main()
