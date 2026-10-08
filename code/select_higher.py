"""Reproduce the bounded seed/conductor search above degree 32.
Only seeds with fundamental class number <=64 are used; no global claim.
"""
from run_paths import output_path,generated_input
import json
from pathlib import Path
from flint import arb,ctx
from catalogue import prime_factors,kronecker_prime
ROOT=Path(__file__).resolve().parents[1]
ctx.dps=40
cat=json.loads(generated_input('catalogue64.json').read_text())
fields=[o for o in cat['orders'] if o['conductor']==1 and o['d']%32==20]
records=[]
for budget in (32,64,80,96,128):
    best=None
    for o in fields:
        dk,hk=o['d'],o['h'];costmax=2*budget//hk
        for f in range(1,2*costmax**2+1,2):
            cost=f
            for p in prime_factors(f):cost=cost//p*(p-kronecker_prime(-dk,p))
            degree=hk*cost//2
            if degree<=budget:
                N=dk*f*f//4
                if best is None or N>best['N']:
                    best=dict(budget=budget,N=N,d=4*N,h=hk*cost,predicted_degree=degree,
                              fundamental_d=dk,conductor=f,
                              rate_estimate=str(arb.pi()*arb(N).sqrt()/arb(10).log()-8*arb(2).log()/arb(10).log()))
    records.append(best)
output_path('higher_candidates_unverified.json').write_text(json.dumps(
    dict(status='candidate selection only; incomplete seed fields and unverified field degrees',rows=records),indent=2))
print([(r['budget'],r['N'],r['predicted_degree']) for r in records])
