"""Pi-free algebraic construction of A,B,x, using a 17-isogeny.

The input A0,B0,x0 is the published Borwein N=253 identity. Phi_17 is
Sutherland's classical modular polynomial. The derivation is in REPORT.md.
Only the final comparison calls arb.pi(); the construction does not.
"""
import json,re,math,time,hashlib
from pathlib import Path
from flint import arb,ctx,fmpz_poly
from run_paths import output_path
ROOT=Path(__file__).resolve().parents[1]

def read_phi():
    out=[]
    for line in (ROOT/'data'/'phi_j_17.txt').read_text().splitlines():
        m=re.fullmatch(r'\[(\d+),(\d+)\]\s+(-?\d+)',line.strip())
        assert m,line
        i,j,c=map(int,m.groups());out.append((i,j,c))
        if i!=j:out.append((j,i,c))
    return out

def derivatives(phi,J,K):
    jp=[J**i for i in range(19)];kp=[K**i for i in range(19)]
    values={}
    for a,b in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2)):
        v=arb(0)
        for i,j,c in phi:
            if i<a or j<b:continue
            c*=math.prod(range(i-a+1,i+1))*math.prod(range(j-b+1,j+1))
            v+=c*jp[i-a]*kp[j-b]
        values[a,b]=v
    return values

def C(t):return (t-512)*(t+64)/((t+256)*(t-64))
def Lrho(t):return 96*t/((t+256)*(t-64))
def j_of_t(t):return (t+256)**3/t**2

def isolated_root(poly):
    """Select the real embedding from FLINT's certified root isolation."""
    hits=[z.real for z,e in poly.complex_roots()
          if e==1 and z.imag.is_zero() and -arb('4e-367')<z.real<-arb('2e-367')]
    if len(hits)!=1:
        raise ArithmeticError('Expected one simple real root in the defining interval')
    return hits[0]


def evaluate(coefficients,x):
    value=arb(0)
    for coefficient in reversed(coefficients):
        value=value*x+coefficient
    return value


def refined_root(poly):
    """Isolate once at 80 digits, then retain a certified root enclosure.

    If r is in I and m is its midpoint, the mean value theorem gives
    r in m - p(m)/p'(I), provided 0 is not in p'(I). Intersecting this
    enclosure with I retains r. Arb evaluates every operation outwardly;
    ordinary floating-point Newton iterates would not suffice here.
    """
    target=ctx.prec
    try:
        ctx.dps=80
        root=isolated_root(poly)
        coefficients=list(poly.coeffs())
        derivative=[i*c for i,c in enumerate(coefficients)][1:]
        while True:
            ctx.prec=min(target,2*ctx.prec)
            for _ in range(8):
                midpoint=root.mid()
                slope=evaluate(derivative,root)
                if slope.contains(0):
                    raise ArithmeticError('Root refinement derivative contains zero')
                root=root.intersection(midpoint-evaluate(coefficients,midpoint)/slope)
                if root.rel_accuracy_bits()>=ctx.prec-12:
                    break
            else:
                raise ArithmeticError('Root refinement did not reach the requested accuracy')
            if ctx.prec==target:
                return root
    finally:
        ctx.prec=target


def construct(dps=2000,root_method='newton'):
    ctx.dps=dps;p=17;N=73117
    Y=2216752650+668376072*arb(11).sqrt();x0=-1/Y**2
    A0=(-37515813+11937508*arb(11).sqrt())/6523272
    B0=-arb(9686105)/543606+arb(8291270)*arb(11).sqrt()/815409
    assert (B0*B0-253*(1-x0)).contains(0)
    s0=B0/arb(253).sqrt()
    # Rationalize t=64(1+s)/(1-s) to avoid subtracting nearly equal numbers.
    t0=64*(1+s0)**2/x0
    J=j_of_t(t0)
    psi0=1+6*(-A0/B0+Lrho(t0))/C(t0)
    data=json.loads((ROOT/'results'/'candidate32.json').read_text())
    poly=fmpz_poly(data['p'])
    if root_method=='newton':x=refined_root(poly)
    elif root_method=='all-roots':x=isolated_root(poly)
    else:raise ValueError(root_method)
    B=(N*(1-x)).sqrt();s=B/arb(N).sqrt()
    t=64*(1+s)**2/x;K=j_of_t(t)
    dv=derivatives(read_phi(),J,K)
    assert dv[0,0].contains(0)
    Pj,Pk,Pjj,Pjk,Pkk=(dv[z] for z in ((1,0),(0,1),(2,0),(1,1),(0,2)))
    assert not Pj.contains(0) and not Pk.contains(0)
    R=-Pj*J/(p*Pk*K)
    dJlogR=Pjj/Pj+1/J-Pjk/Pk
    dKlogR=Pjk/Pj-Pkk/Pk-1/K
    L=-J*dJlogR-p*R*K*dKlogR
    f0=psi0/6-J/(2*(J-1728))+arb(1)/3
    psi1=6*(L+f0)/(p*R)+3*K/(K-1728)-2
    A=B*((1-psi1)*C(t)/6+Lrho(t))
    return A,B,x,dict(J=str(J),K=str(K),R=str(R),psi0=str(psi0),psi1=str(psi1))

def main():
    start=time.perf_counter();A,B,x,extra=construct()
    assert A>0 and B>0 and -1<x<0
    # c_n=(1/4)_n(1/2)_n(3/4)_n/(n!)^3, each recurrence factor <1.
    # For n>=m the absolute term ratio is bounded by |x|*(1+B/(A+Bm)).
    results=[];c=arb(1);power=arb(1);total=arb(0)
    for n in range(4):
        total+=(A+B*n)*c*power
        m=n+1
        c*=((arb(n)+arb(1)/4)*(arb(n)+arb(1)/2)*(arb(n)+arb(3)/4))/(n+1)**3
        power*=x
        first=abs((A+B*m)*c*power)
        ratio=abs(x)*(1+B/(A+B*m));assert ratio<1
        tail=first/(1-ratio)
        enclosure=total+arb(0,tail.upper())
        pi_interval=1/enclosure
        # An independent numerical comparison, not used to choose the coefficients.
        assert pi_interval.contains(arb.pi())
        approximation=(1/total).mid()
        error=abs(pi_interval.mid()-approximation)+pi_interval.rad()
        error_digits=int(float(-error.upper().log()/arb(10).log()))
        assert error<arb(10)**(-error_digits)
        results.append(dict(terms=m,tail_1_over_pi=str(tail),pi_radius=str(pi_interval.rad()),
                            approximation_error_bound=str(error),certified_absolute_error_digits=error_digits))
    result=dict(N=73117,rate='366.5212453310619038906676419874',
                pi_used_to_construct_coefficients=False,phi17_source='https://math.mit.edu/~drew/modpolys/jfiles/phi_j_17.txt',
                phi17_sha256=hashlib.sha256((ROOT/'data'/'phi_j_17.txt').read_bytes()).hexdigest(),
                A=str(A),B=str(B),x=str(x),truncations=results,elapsed=time.perf_counter()-start)
    output_path('identity32.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(N=result['N'],pi_used_to_construct_coefficients=False,
                         truncations=[{'terms':r['terms'],'error_digits':r['certified_absolute_error_digits']} for r in results],
                         elapsed=result['elapsed']),indent=2))

if __name__=='__main__':main()
