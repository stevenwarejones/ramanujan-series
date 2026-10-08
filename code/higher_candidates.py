"""Exact parameter/radical checks, NOT optimality certificates, above degree 32."""
from run_paths import output_path,generated_input
import json,time
from pathlib import Path
from flint import fmpz_poly,arb,acb,ctx
from extend32 import low_degree_cover,remove_zero
from new_candidate import radical_certificate
from polynomials import canon
ROOT=Path(__file__).resolve().parents[1]

def main():
    start=time.perf_counter()
    rows=json.loads(generated_input('higher_candidates_unverified.json').read_text())['rows']
    out=[]
    for row in rows:
        if row['budget']==32:continue
        N=row['N'];d=row['d'];m=row['predicted_degree'];ctx.dps=160
        H=fmpz_poly.hilbert_class_poly(-d);assert H.degree()==row['h']
        q=remove_zero(low_degree_cover(H,2,m));content,fs=q.factor();check=fmpz_poly([content])
        for p,e in fs:check*=p**e
        assert q==check
        tau=acb(arb(-1)/2,arb(N).sqrt()/2)
        t=(tau.modular_eta()/(2*tau).modular_eta())**24;x=256*t/(t+64)**2
        matches=[]
        for p,e in fs:
            val=acb(0)
            for c in reversed(p.coeffs()):val=val*x+int(c)
            if val.contains(0):matches.append(canon(p))
        assert len(matches)==1
        p=matches[0];assert p.degree()==m
        cert=radical_certificate(p,N);assert cert['full_x_B_degree']==m
        rate=-abs(x).log()/arb(10).log()
        r=dict(row,degree=m,rate=str(rate),p=list(map(int,p.coeffs())),radical=cert,
               optimality='unproved; seed search used only fundamental fields of class number <=64',
               constant_A='exists in Q(x) by the CM descent proof; not expanded here',
               parameter_and_B_status='exact irreducible polynomial and square identity checked')
        out.append(r)
        output_path('higher_candidates.json').write_text(json.dumps({'rows':out},indent=2))
        print('HIGHER',row['budget'],N,m,'rate',str(rate)[:70],'seconds',round(time.perf_counter()-start,1),flush=True)

if __name__=='__main__':main()
