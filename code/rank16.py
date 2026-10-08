from run_paths import output_path
import json,hashlib,time
from pathlib import Path
from flint import fmpz_poly,arb,ctx
OUT=Path(__file__).resolve().parents[1]/'data'
r=json.loads((OUT/'search_results.json').read_text());ctx.dps=120
start=time.perf_counter();cache={};exclusions=0;rootsets=0

def roots(p):
 global rootsets
 key=tuple(p)
 if key not in cache:
  z=fmpz_poly(p).complex_roots();assert sum(e for _,e in z)==len(p)-1
  assert all(e==1 for _,e in z)
  cache[key]=[w for w,e in z];rootsets+=1
 return cache[key]

wins={}
for N in (2,3,4):
 for k in range(1,17):
  rec=r['frontiers'][str(N)][k];zs=roots(rec['p']);small=None
  for z in zs:
   if all(z is zz or abs(z)<abs(zz) for zz in zs):small=z;break
  assert small is not None,(N,k)
  wins[N,k]=(rec,small)

for budget in range(1,17):
 for N in (3,4):
  assert abs(wins[2,budget][1]) < abs(wins[N,budget][1])

for rec in r['records']:
 N=rec['level'];bd=rec['B_degree'];p=rec['p'];m=len(p)-1
 wr,wz=wins[N,bd]
 # Exact integer Rouche exclusion at a disk strictly larger than the winner.
 R=-abs(wz).log()/arb(10).log();k=int(float(R))
 assert arb(k)<R and R<arb(k+1)
 den=10**k
 assert abs(wz)<arb(1)/den
 sm=0
 for coef in reversed(p[1:]):sm=sm*den+abs(coef)
 # sum |p_i| den^(m-i), i=1..m, calculated directly (reverse horner above differs)
 sm=sum(abs(p[i])*den**(m-i) for i in range(1,m+1))
 if abs(p[0])*den**m>sm:
  exclusions+=1;continue
 zz=roots(p)
 for budget in range(bd,17):
  old,ow=wins[N,budget]
  if p==old['p']:continue
  assert all(abs(ow)<abs(z) for z in zz),(N,bd,budget,rec['d'],old['d'])

res={'candidate_polynomials':len(r['records']),'exact_Rouche_exclusions':exclusions,'ball_root_sets':rootsets,
 'budgets_certified':list(range(1,17)),'elapsed_seconds':time.perf_counter()-start,
 'sha256_search_results':hashlib.sha256((OUT/'search_results.json').read_bytes()).hexdigest(),
 'winners':{str(N):{'d':wins[N,16][0]['d'],'parameter_degree':wins[N,16][0]['parameter_degree'],
 'full_degree':wins[N,16][0]['B_degree'],'x':str(wins[N,16][1]),
 'rate':str(-abs(wins[N,16][1]).log()/arb(10).log())} for N in (2,3,4)}}
output_path('reranked16.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
