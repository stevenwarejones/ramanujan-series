#!/usr/bin/env python3
"""Exact finite functional reduction test; says nothing about CM-value independence.

The target is C(q)=-D^-2[(D^-1 V(q))^2], V=16 A_3,
A_r=sum (sigma_r(n)-sigma_r(n/2))q^n. Test the Q-linear span of
1, D^-j A_r(q^d) for r=3,5,7; d=1,2; j=1,2,3,4;
and the three products R(q)^2, R(q)R(q^2), R(q^2)^2, R=D^-2 V.
A rational separating functional on q^0,...,q^60 proves nonmembership.
"""
import json
from fractions import Fraction as F
from math import comb


def rref(rows):
    a=[list(map(F,row)) for row in rows]; rank=0; pivots=[]
    for col in range(len(a[0])-1):
        pivot=next((i for i in range(rank,len(a)) if a[i][col]),None)
        if pivot is None: continue
        a[rank],a[pivot]=a[pivot],a[rank]
        d=a[rank][col]; a[rank]=[v/d for v in a[rank]]
        for i in range(len(a)):
            if i!=rank and a[i][col]:
                d=a[i][col];a[i]=[v-d*w for v,w in zip(a[i],a[rank])]
        pivots.append(col);rank+=1
        if rank==len(a): break
    return a,pivots


def build(cutoff=60):
    kernels={}
    for power in (3,5,7):
        sigma=[0]*(cutoff+1)
        for d in range(1,cutoff+1):
            for n in range(d,cutoff+1,d):sigma[n]+=d**power
        kernels[power]=[0]+[sigma[n]-(sigma[n//2] if n%2==0 else 0) for n in range(1,cutoff+1)]
    labels=['1']; columns=[[F(1)]+[F(0)]*cutoff]
    for r in (3,5,7):
        for d in (1,2):
            for j in (1,2,3,4):
                labels.append(f'D^-{j} A_{r}(q^{d})')
                columns.append([F(0)]+[F(kernels[r][n//d],n**j) if n%d==0 else F(0) for n in range(1,cutoff+1)])
    rs={d:[F(0)]+[F(16*kernels[3][n//d],(n//d)**2) if n%d==0 else F(0) for n in range(1,cutoff+1)] for d in (1,2)}
    for d,e in ((1,1),(1,2),(2,2)):
        labels.append(f'R(q^{d}) R(q^{e})')
        columns.append([sum((rs[d][k]*rs[e][n-k] for k in range(n+1)),F(0)) for n in range(cutoff+1)])
    p=[F(0)]+[F(16*kernels[3][n],n) for n in range(1,cutoff+1)]
    target=[F(0)]+[-sum((p[k]*p[n-k] for k in range(1,n)),F(0))/n**2 for n in range(1,cutoff+1)]
    for n in range(1,cutoff+1):
        direct=sum((F(256*comb(n//k+1,3),k**4) for k in range(1,n+1,2) if n%k==0 and n//k>=2),F(0))
        reduced=F(128,3*n**4)*(kernels[7][n]-kernels[5][n])
        assert direct==reduced
    return labels,columns,target


def main():
    labels,columns,target=build()
    # Solve for a linear functional killing the proposed basis and sending C to 1.
    reduced,pivots=rref([col+[F(0)] for col in columns]+[target+[F(1)]])
    assert all(any(row[:-1]) or not row[-1] for row in reduced)
    weights=[F(0)]*len(target)
    for row,pivot in zip(reduced,pivots):weights[pivot]=row[-1]
    assert all(sum((a*b for a,b in zip(weights,col)),F(0))==0 for col in columns)
    assert sum((a*b for a,b in zip(weights,target)),F(0))==1
    print(json.dumps({'target':'C(q)=-D^-2[(D^-1 V)^2]', 'maximum_q_degree':60,
      'basis':labels,'basis_rank':len(pivots)-1,'augmented_rank':len(pivots),
      'witness':{n:str(v) for n,v in enumerate(weights) if v},
      'witness_on_each_basis_element':'0','witness_on_target':'1',
      'independent_phase_reduction':'sum_(k odd) a_k(q)^4=(128/3) D^-4(A_7-A_5)',
      'scope':'No Q-linear functional identity within this stated 28-function ansatz. Does NOT exclude other functions, algebraic coefficients, nonlinear relations, or a relation only at q=exp(-pi).'},indent=2))

if __name__=='__main__':main()
