"""Exact q-series/Sturm certificate for the external Phi_17 data.

Delta(tau)^18 Delta(17tau)^18 Phi(j(tau),j(17tau)) has weight 432
on Gamma_0(17), index 18, hence Sturm bound 648. The Delta prefactor
starts with q^324. F=q^306 Phi is holomorphic at infinity, so checking
F through q^630 suffices. All arithmetic below is over the integers.
"""
import json,time
from pathlib import Path
from flint import fmpz_poly,fmpz_series,ctx
from identity32 import read_phi
ROOT=Path(__file__).resolve().parents[1]

def main():
    start=time.perf_counter();precision=631;ctx.cap=precision+2
    a=[0]*(precision+2);b=[0]*(precision+2);a[0]=b[0]=1
    for d in range(1,precision+2):
        for k in range(d,precision+2,d):a[k]+=240*d**3;b[k]-=504*d**5
    e4=fmpz_series(a);e6=fmpz_series(b)
    diff=e4**3-e6**2
    delta_over_q=[int(diff[i+1])//1728 for i in range(precision)]
    assert all(int(diff[i+1])%1728==0 for i in range(precision))
    assert delta_over_q[0]==1
    J=(e4**3)/fmpz_series(delta_over_q,prec=precision)
    jpoly=fmpz_poly([int(J[i]) for i in range(precision)])
    jk=[0]*precision
    for i in range((precision-1)//17+1):jk[17*i]=int(J[i])
    kp=fmpz_poly(jk)
    jpows=[fmpz_poly([1])];kpows=[fmpz_poly([1])]
    for k in range(18):
        jpows.append((jpows[-1]*jpoly).truncate(precision))
        kpows.append((kpows[-1]*kp).truncate(precision))
    result=fmpz_poly([])
    for i,j,c in read_phi():
        shift=306-i-17*j;assert shift>=0
        term=(jpows[i]*kpows[j]).truncate(precision-shift)
        result+=(fmpz_poly([0]*shift+list(term.coeffs()))*c)
    assert not result
    out=dict(status='PASS',integer_q_coefficients_checked=631,weight=432,
             group='Gamma_0(17)',group_index=18,Sturm_bound=648,
             elapsed=time.perf_counter()-start)
    (ROOT/'results'/'phi17_check.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out))

if __name__=='__main__':main()
