"""Independent finite catalogue; no recursive conductor search.

Completeness beyond the finite enumeration relies on Watkins (2004), p.23.
For f>1, phi(f) <= f prod(1-chi(p)/p) <= u H/h_K, and
phi(f) >= sqrt(f/2), so f <= 2 (u H/h_K)^2 suffices.
"""
from run_paths import output_path
import argparse, json, math, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WATKINS = [
 (9,163),(18,427),(16,907),(54,1555),(25,2683),(51,3763),(31,5923),(131,6307),
 (34,10627),(87,13843),(41,15667),(206,17803),(37,20563),(95,30067),(68,34483),(322,31243),
 (45,37123),(150,48427),(47,38707),(350,58507),(85,61483),(139,85507),(68,90787),(511,111763),
 (95,93307),(190,103027),(93,103387),(457,126043),(83,166147),(255,134467),(73,133387),(708,164803),
 (101,222643),(219,189883),(103,210907),(668,217627),(85,158923),(237,289963),(115,253507),(912,260947),
 (109,296587),(339,280267),(106,300787),(691,319867),(154,308323),(268,462883),(107,375523),(1365,335203),
 (132,393187),(345,389467),(159,546067),(770,439147),(114,425107),(427,532123),(163,452083),(1205,494323),
 (179,615883),(291,586987),(128,474307),(1302,662803),(132,606643),(323,647707),(216,991027),(1672,693067)]

def kronecker_prime(D,p):
    if D % p == 0: return 0
    if p == 2: return 1 if D % 8 in (1,7) else -1
    return 1 if pow(D % p, (p-1)//2, p)==1 else -1

def prime_factors(n):
    result=[];p=2
    while p*p<=n:
        if n%p==0:
            result.append(p)
            while n%p==0:n//=p
        p=3 if p==2 else p+2
    if n>1:result.append(n)
    return result

def reduced_forms(d):
    """All primitive reduced forms with discriminant -d, signed b."""
    out=[]
    for a in range(1,math.isqrt(d//3)+1):
        for b in range(-a,a+1):
            if (b*b+d)%(4*a):continue
            c=(b*b+d)//(4*a)
            if c<a or math.gcd(a,math.gcd(b,c))!=1:continue
            if (abs(b)==a or a==c) and b<0:continue
            out.append((a,b,c))
    return out

def build(H):
    assert 1<=H<=len(WATKINS)
    t=time.perf_counter();limit=max(x[1] for x in WATKINS[:H])
    sf=bytearray(b'\x01')*(limit+1)
    for p in range(2,math.isqrt(limit)+1):
        if sf[p]:sf[p*p::p*p]=b'\x00'*(limit//(p*p))
    fund=bytearray(limit+1)
    for d in range(3,limit+1):
        fund[d]=(d%4==3 and sf[d]) or (d%4==0 and sf[d//4] and d//4%4 in (1,2))
    counts=[0]*(limit+1)
    # Sweep a,b and c, independently of the later conductor enumeration.
    for a in range(1,math.isqrt(limit//3)+1):
        for b in range(a+1):
            for c in range(a,(limit+b*b)//(4*a)+1):
                d=4*a*c-b*b
                if fund[d]:counts[d]+=1 if b==0 or b==a or a==c else 2
    fields=[(d,counts[d]) for d in range(3,limit+1) if fund[d] and counts[d]<=H]
    for h,(num,mx) in enumerate(WATKINS[:H],1):
        ds=[d for d,hh in fields if hh==h]
        assert len(ds)==num and max(ds)==mx,(h,len(ds),max(ds))
    maxf=2*(3*H)**2
    factors=[prime_factors(f) for f in range(maxf+1)]
    orders={}
    for dk,hk in fields:
        units=3 if dk==3 else 2 if dk==4 else 1
        for f in range(1,2*(units*H//hk)**2+1):
            cost=f
            for p in factors[f]:cost=cost//p*(p-kronecker_prime(-dk,p))
            if f==1:h=hk
            else:
                assert hk*cost%units==0
                h=hk*cost//units
            if h<=H:
                d=dk*f*f
                assert d not in orders
                orders[d]=dict(d=d,h=h,fundamental_d=dk,conductor=f)
    result=dict(class_number_bound=H,discriminant_bound=limit,fields=len(fields),
                orders=[orders[d] for d in sorted(orders)],watkins_checks=H,
                conductor_bound='f <= 2 (u H/h_K)^2',elapsed=time.perf_counter()-t)
    path=output_path(f'catalogue{H}.json');path.write_text(json.dumps(result,indent=2))
    print('CATALOGUE',H,len(fields),len(orders),'max order',max(orders),'seconds',result['elapsed'],flush=True)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--bound',type=int,default=32)
    args=ap.parse_args();r=build(args.bound)
    if args.bound==32:
        old=json.loads((ROOT/'data'/'orders32.json').read_text())
        assert {x['d']:(x['h'],x['fundamental_d'],x['conductor']) for x in r['orders']} == {x['d']:(x['h'],x['fundamental_d'],x['conductor']) for x in old['orders']}
        print('Independent conductor enumeration matches every legacy order.')
