#!/usr/bin/env python3
"""Certified harmonic-balance evaluation on the square torus; stdlib only.

Proof and tail estimates are in NOTE.md. No output files
are changed by running this script. Fourier integration is performed by
exact frequency matching with outward interval coefficients, not a grid.
"""
import json
from math import factorial
from intervals import I, SCALE, decimal_string, pi_interval, exp_positive, gauss_log_jet, log_jet


def conv(left,right):
    zero=I.exact(0); result={}
    for i,a in left.items():
        for j,b in right.items():
            result[i+j]=result.get(i+j,zero)+a*b
    return result


def moment_even(poly):
    return sum((a*poly[-n] for n,a in poly.items()),I.exact(0))


def finite_moments(a,degree=6):
    wave={n:v for k,v in a.items() for n in (k,-k)}
    square=conv(wave,wave)
    moments={2:moment_even(wave),4:moment_even(square)}
    if degree>=6:
        cube=conv(square,wave)
        moments[6]=moment_even(cube)
    return moments


def enclose(v,bound):
    return I(v.lo-bound.hi,v.hi+bound.hi)


def encode(v):
    return {'lower':decimal_string(v.lo),'upper':decimal_string(v.hi),
            'width':decimal_string(v.hi-v.lo)}


def main():
    q=exp_positive(pi_interval()).reciprocal(); p=q.sqrt()
    cutoff=241
    a={k:4*p**k/(k*(1-q**k)) for k in range(1,cutoff+1,2)}
    moments=finite_moments(a)
    N=cutoff+2
    epsilon=8*p**N/(N*(1-q**N)*(1-q))
    size=2*sum(a.values(),I.exact(0))
    raw=[I.exact(1)]+[I.exact(0)]*6
    for j,m in moments.items():
        error=j*(size+epsilon)**(j-1)*epsilon
        raw[j]=enclose(m,error)/factorial(j)
    logs=log_jet(raw)
    reference=gauss_log_jet(terms=650,order=6)
    for j in (2,4,6):
        assert max(logs[j].lo,reference[j].lo)<=min(logs[j].hi,reference[j].hi)
        assert logs[j].hi-logs[j].lo<SCALE//10**150
    approx={}
    for count in (1,2,3,5,10):
        subset={k:a[k] for k in range(1,2*count,2)}
        mm=finite_moments(subset,4)
        cc=(mm[4]-3*mm[2]**2)/24
        approx[str(count)]={'highest_frequency':2*count-1,'c':encode(cc),
                            'true_minus_truncated':encode(logs[4]-cc)}
    first=-a[1]**4/4
    two=-(a[1]**4+a[3]**4)/4+a[1]**3*a[3]/3
    three=-(a[1]**4+a[3]**4+a[5]**4)/4+a[1]**3*a[3]/3+a[1]**2*a[3]*a[5]+a[1]*a[3]**2*a[5]
    for target,count in [(first,1),(two,2),(three,3)]:
        mm=finite_moments({k:a[k] for k in range(1,2*count,2)},4)
        residual=target-(mm[4]-3*mm[2]**2)/24
        assert residual.lo<=0<=residual.hi
    # Exact all-mode independent-phase coefficient, with a positive tail.
    independent=-sum((v**4 for v in a.values()),I.exact(0))/4
    independent_tail=64*q**(2*N)/(N**4*(1-q**N)**4*(1-q**4))
    independent=I(independent.lo-independent_tail.hi,independent.hi)
    resonance=a[1]**3*a[3]/3
    correction=logs[4]-first
    fraction=resonance/correction
    assert correction.lo>0 and fraction.lo>0 and fraction.hi<SCALE
    # Frequency balance improves the quartic tail from p^N to q^N.
    # Monotonicity of consecutive odd-harmonic truncations makes it one-sided.
    assert q.hi<SCALE//20
    fast_cutoff=121; NN=fast_cutoff+2
    small=finite_moments({k:v for k,v in a.items() if k<=fast_cutoff},4)
    small_c=(small[4]-3*small[2]**2)/24
    C=4/(1-q)
    bracket=(NN+2)**2/(1-q)+2*(NN+2)*q/(1-q)**2+q*(1+q)/(1-q)**3
    fast_tail=2*C**4*q**NN*bracket/3
    fast_c=I(small_c.lo,small_c.hi+fast_tail.hi)
    assert max(fast_c.lo,reference[4].lo)<=min(fast_c.hi,reference[4].hi)
    assert fast_c.hi-fast_c.lo<SCALE//10**150
    print(json.dumps({'method':'integer outward intervals; exact frequency matching; uniform Fourier tail',
      'positive_odd_harmonics':len(a),'highest_frequency':cutoff,
      'uniform_log_wave_tail':encode(epsilon),
      'log_coefficients':{j:encode(logs[j]) for j in (2,4,6)},
      'independent_gauss_intersections':True,'all_widths_below':'1e-150',
      'truncations':approx,'independent_phase_c':encode(independent),
      'leading_1_plus_1_plus_1_equals_3_resonance':encode(resonance),
      'fraction_of_single_harmonic_correction':encode(fraction),
      'frequency_balance_certificate':{'positive_odd_harmonics':61,
          'highest_frequency':fast_cutoff,'first_omitted_frequency':NN,
          'one_sided_tail_bound':encode(fast_tail),'c':encode(fast_c),
          'independent_gauss_intersection':True,'width_below':'1e-150'},
      'scope':'Certified evaluation and decomposition, not a familiar-constant reduction or novelty proof.'},indent=2))

if __name__=='__main__': main()
