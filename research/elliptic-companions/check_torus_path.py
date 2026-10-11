#!/usr/bin/env python3
"""Optional numerical cross-check of the lifted path and theta product.

Requires mpmath. These are floating checks, separate from verify_fourier.py's
rigorous interval certificates. Analytic derivations are in NOTE.md.
"""
import json
import mpmath as mp
mp.mp.dps=100
r=1+mp.sqrt(2); alpha=r*r; beta=1/alpha
K=mp.ellipk(mp.mpf('.5')); omega=K/mp.sqrt(2); q=mp.exp(-mp.pi)
errors=[]; rows=[]
for numerator in (0,1,2,3,4):
    x=r*numerator/4
    y=mp.sqrt((1+x*x)*(alpha+beta*x*x))
    t=mp.quad(lambda xx:1/mp.sqrt((1+xx*xx)*(alpha+beta*xx*xx)),[0,x])
    z=1j*omega/2-t
    sn=mp.ellipfun('sn',mp.sqrt(2)*z,mp.mpf('.5'))
    cn=mp.ellipfun('cn',mp.sqrt(2)*z,mp.mpf('.5'))
    dn=mp.ellipfun('dn',mp.sqrt(2)*z,mp.mpf('.5'))
    XX=-1+2/sn**2; YY=-2*mp.sqrt(2)*cn*dn/sn**3
    X=mp.sqrt(2)*(x+1j)/(x-1j)-1; Y=mp.sqrt(2)*1j*y/(x-1j)**2
    theta=mp.pi*t/omega
    logT=mp.log((alpha+beta*x*x)/(1+x*x))
    theta_log=2*mp.log(mp.jtheta(3,theta/2,mp.sqrt(q))/mp.jtheta(4,theta/2,mp.sqrt(q)))
    fourier=8*mp.fsum(q**(mp.mpf(k)/2)*mp.cos(k*theta)/(k*(1-q**k)) for k in range(1,161,2))
    error=max(abs(XX-X),abs(YY-Y),abs(logT-theta_log),abs(logT-fourier))
    assert error<mp.mpf('1e-94')
    errors.append(error)
    rows.append({'x_over_r':str(mp.mpf(numerator)/4),'theta_over_pi':mp.nstr(theta/mp.pi,80),'maximum_error':mp.nstr(error,6)})
assert abs(t-omega/2)<mp.mpf('1e-95')
# Independent all-order integral identity at noninfinitesimal parameters.
checks=[]
for s in (mp.mpf('-0.7'),mp.mpf('0.3'),mp.mpf('1.2')):
    wave=lambda th:2*mp.log(mp.jtheta(3,th/2,mp.sqrt(q))/mp.jtheta(4,th/2,mp.sqrt(q)))
    average=mp.quad(lambda th:mp.cosh(s*wave(th)),[0,mp.pi/4,mp.pi/2])*2/mp.pi
    original=2**(-s)*mp.hyp2f1(mp.mpf('.5')+s,mp.mpf('.5')+s,1,mp.mpf('.5'))/mp.hyp2f1(mp.mpf('.5'),mp.mpf('.5'),1,mp.mpf('.5'))
    error=abs(average-original)
    assert error<mp.mpf('1e-94')
    checks.append({'s':str(s),'absolute_error':mp.nstr(error,6)})
# Canonical translated first-Kronecker alphabet, including its constant.
u0=(-1+1j)/4
bs=[(3-1j)/4,(3+1j)/4,(1+1j)/4,(1-1j)/4]
alphabet_errors=[]
for v in (mp.mpf(0),mp.mpf(1)/16,mp.mpf(1)/8,mp.mpf(1)/4):
    u=u0+v
    eta=mp.pi*sum(sign*mp.jtheta(j,mp.pi*u,q,1)/mp.jtheta(j,mp.pi*u,q)
                  for sign,j in ((1,2),(1,3),(-1,1),(-1,4)))
    mapped=sum(sign*mp.pi*mp.jtheta(1,mp.pi*(v-b),q,1)/mp.jtheta(1,mp.pi*(v-b),q)
               for sign,b in zip((1,1,-1,-1),bs))
    error=abs(eta-mapped)
    assert error<mp.mpf('1e-94')
    alphabet_errors.append(mp.nstr(error,6))
print(json.dumps({'precision':mp.mp.dps,'canonical_alphabet_errors':alphabet_errors,'scope':'Floating-point cross-checks, not interval certificates.','path_checks':rows,'generating_function_checks':checks},indent=2))
