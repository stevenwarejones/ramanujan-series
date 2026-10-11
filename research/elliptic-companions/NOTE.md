# An elliptic logarithmic coefficient and its Fourier interactions

This note gives explicit proofs and reproducible certificates for one
hypergeometric coefficient. It makes no novelty or irrationality claim.
The [bounded reduction audit](BOUNDED_REDUCTION.md) records a separate,
limited attempt to express it through simpler quantities.

## 1. Target and normalization

Write F for the Gauss hypergeometric function and define

```math
h(s)=2^{-s}\frac{F(\tfrac12+s,\tfrac12+s;1;\tfrac12)}
{F(\tfrac12,\tfrac12;1;\tfrac12)},\qquad
\log h(s)=R s^{2}+c s^{4}+\ell_3 s^{6}+\cdots.
```

The affine normalization does not change coefficients of degree at least two.
All coefficients are ordinary Taylor coefficients, not derivatives. Elliptic
integrals use the **parameter** convention

```math
K(m)=\int_0^{\pi/2}\frac{d\phi}{\sqrt{1-m\sin^{2}\phi}},
\qquad K=K(1/2)=\frac{\Gamma(1/4)^{2}}{4\sqrt{\pi}}.
```

The classical Legendre representation and Pfaff transformation [1] give

```math
h(s)=\frac{\int_0^\pi A(\phi)^{-1/2}e^{s\log A(\phi)}\,d\phi}
{\int_0^\pi A(\phi)^{-1/2}\,d\phi},\qquad
A(\phi)=3+2\sqrt2\cos\phi.
```

One can verify the normalization from
P_{−1/2+s}(3)=2^{−1/2−s}F(1/2+s,1/2+s;1;1/2).
Put r=1+√2, α=r², β=1/α and set x=tan(φ/2). Then

```math
y^{2}=(1+x^{2})(\alpha+\beta x^{2}),\qquad
T(x)=\frac{\alpha+\beta x^{2}}{1+x^{2}},\qquad
\frac{d\phi}{\sqrt{A(\phi)}}=\frac{2\,dx}{y}.
```

Take the positive y on the real path. The involution x↦α/x sends T↦1/T
and dx/y↦−dx/y, with fixed point x=r. Odd moments vanish. For even j,

```math
m_j=4\int_0^{r}\log(T(x))^{j}\frac{dx}{y},\qquad m_0=\sqrt2K,
\qquad R=\frac{m_2}{2m_0},\qquad
c=\frac{m_4}{24m_0}-\frac{m_2^{2}}{8m_0^{2}}.
```

The value of m₀ follows from the classical quadratic transformation for K [2].

## 2. The square elliptic curve and the actual path

The exact map

```math
X=\sqrt2\frac{x+i}{x-i}-1,\qquad
Y=\frac{\sqrt2\,i\,y}{(x-i)^{2}}
```

satisfies

```math
E:\ Y^{2}=X^{3}-X,\qquad \frac{dX}{Y}=-2\frac{dx}{y},\qquad
 g=-\frac{iY}{X+1}=\frac{y}{1+x^{2}},\qquad T=g^{2}.
```

All identities are algebraic substitutions. The logarithmic function has divisor
(g)=(0,0)+(1,0)−(−1,0)−O, so its zeros and poles are points of order two.
The endpoints are

```math
P_0=(-r,-i\sqrt2r),\quad 2P_0=(-1,0),\qquad
P_1=(i,-1+i),\quad 2P_1=(0,0).
```

The doubling formula proves both endpoints have exactly order four.
`audit_geometry.py` checks the curve, differential, function, involution,
endpoints, and doublings symbolically.

Set ω=K/√2. In the Weierstrass uniformization X=℘(z;4,0), Y=℘′(z;4,0)/2,
the lattice is 2ω(ℤ+iℤ). The lift is fixed, including its winding:

```math
z(x)=\frac{i\omega}{2}-\int_0^{x}\frac{dt}{y(t)}.
```

Indeed ℘(ω/2)=r by duplication, while ℘(iz)=−℘(z) and ℘′(iz)=i℘′(z)
fix the branch at P₀. The differential fixes the rest of the lift, and
∫₀ʳ dx/y=m₀/4=ω/2. Thus the path γ from P₁ to P₀ is the straight segment
from (−1+i)/4 to i/4 in coordinates z/(2ω). Torsion endpoints alone would
not have specified this path.

With ω₀=dx/y and η=dlog g, define the word convention

```math
I_j=\int_\gamma\omega_0(x)\,
\frac{\left(\int_{P_1}^{x}\eta\right)^{j}}{j!}.
```

Here g(P₁)=1 and g is positive along the path, so the inner integral is the
real log g, with no omitted additive constant. Directly,

```math
m_{2k}=-4\,2^{2k}(2k)!I_{2k},\qquad
R=4\frac{I_2}{I_0},\qquad
c=16\frac{I_4}{I_0}-8\left(\frac{I_2}{I_0}\right)^{2}.
```

## 3. The theta quotient and a rule for every coefficient

Let q=e^{−π}, p=√q and θ=(π/ω)∫₀ˣdt/y(t). Standard elliptic/theta
relations and products [2,3] yield on the path

```math
g=\frac{\theta_3(\theta/2,p)}{\theta_4(\theta/2,p)},\qquad
L(\theta)=\log T(x(\theta))
=2\sum_{k\ge1,\ k\ \mathrm{odd}}a_k\cos(k\theta),
\qquad a_k=\frac{4q^{k/2}}{k(1-q^{k})}
=\frac{2}{k\sinh(\pi k/2)}.
```

For completeness, in uniformizing coordinates v=πz/(2ω), the elliptic
relations give g=iθ₂(v,q)θ₃(v,q)/(θ₁(v,q)θ₄(v,q)). Set
E=e^{2iv}=p e^{−iθ}. The factor i cot v=(1+E)/(1−E) and the remaining
products combine into

```math
g=\prod_{n\ge0}
\frac{(1+q^{n+1/2}e^{i\theta})(1+q^{n+1/2}e^{-i\theta})}
{(1-q^{n+1/2}e^{i\theta})(1-q^{n+1/2}e^{-i\theta})}.
```

Each paired factor is positive for real θ. Expanding the real logarithms
proves the displayed Fourier series, absolutely and uniformly. The amplitudes
are positive and decreasing. Also L(−θ)=L(θ), L(θ+π)=−L(θ). Changing variables
in the moments therefore proves, for all complex s,

```math
h(s)=\frac{1}{2\pi}\int_0^{2\pi}e^{sL(\theta)}\,d\theta
=\mathrm{CT}\exp\left(s\sum_{k>0,\ k\ \mathrm{odd}}
 a_k(z^{k}+z^{-k})\right).
```

CT means constant Fourier coefficient on the unit circle. Bounded L and
uniform convergence justify termwise expansion on every bounded s-set.
At degree j, retain exactly the multisets satisfying

```math
\sum_k k(n_k^+-n_k^-)=0,\qquad
\sum_k(n_k^++n_k^-)=j,
```

with weight ∏ₖ aₖ^{nₖ⁺+nₖ⁻}/(nₖ⁺! nₖ⁻!). Parity excludes odd j. At
quadratic order only opposite frequencies survive; at fourth order distinct
frequencies can balance.

For M retained harmonics, the classical Bessel generating function [4] gives

```math
h_M(s)=\sum_{m_1+3m_3+\cdots+(2M-1)m_{2M-1}=0}
\prod_{k=1,3,\ldots,2M-1} I_{m_k}(2a_ks).
```

The limit is locally uniform in s. In particular h₁=I₀(2a₁s), while
h₂=∑ₘ I_{3m}(2a₁s)Iₘ(2a₃s). Independent phases would instead give a product
of I₀ factors, discarding the balanced interactions.

## 4. Quartic decomposition

Write S=∑aₖ², bₙ=∑_{i+j=n}aᵢaⱼ and T₃=∑_{i,j,k}aᵢaⱼaₖa_{i+j+k},
where all indices inside sums are positive odd integers. Then

```math
R=S,\qquad
c=\frac14\sum_n b_n^{2}+\frac13T_3-\frac12S^{2}.
```

Proof: the fourth moment has six sign assignments with two plus and two
minus signs, and eight assignments with three against one. The logarithm
subtracts three times the square of the second moment, divided by 24.
All exchanges of sums are valid by absolute convergence.

Let U sum aᵢaⱼaₖaₗ over **ordered** quadruples with i+j=k+l but unequal
multisets {i,j} and {k,l}. The trivial pairings contribute 2S²−∑aₖ⁴, so

```math
c=-\frac14\sum_k a_k^{4}+\frac14U+\frac13T_3.
```

Both interaction terms are positive. The first term alone is the answer
for independent phases. Keeping the first one, two, or three odd harmonics gives

```math
c_1=-a_1^{4}/4,
```
```math
c_2=-(a_1^{4}+a_3^{4})/4+a_1^{3}a_3/3,
```
```math
c_3=-(a_1^{4}+a_3^{4}+a_5^{4})/4+a_1^{3}a_3/3
+a_1^{2}a_3a_5+a_1a_3^{2}a_5.
```

The last interactions correspond to 1+1+3=5 and 1+5=3+3.
The certified value of (a₁³a₃/3)/(c−c₁) is
0.99891498057706324721107177617… . This is a ratio of corrections.

## 5. Monotonicity and the sharper tail

**Proposition.** For consecutive odd harmonics, c_M increases strictly to c.
For N=2M+1, the first omitted frequency, and 0<q≤1/20, set C=4/(1−q). Then

```math
0<c-c_M\le\frac23C^{4}q^{N}
\left(\frac{(N+2)^{2}}{1-q}
+\frac{2(N+2)q}{(1-q)^{2}}+
\frac{q(1+q)}{(1-q)^{3}}\right).
```

The claim applies to the Fourier-defined coefficient at every such q;
our geometric target has q=e^{−π}.

**Proof of monotonicity.** Adding frequency 3 contributes at least
 a₁³a₃/3−a₃⁴/4>0. Adding odd N≥5 contributes at least
 a₁²a_{N−2}a_N−a_N⁴/4>0, using the three permutations of (1,1,N−2)
in T₃ and decreasing amplitudes. All other new interaction terms are
positive. Absolute convergence gives the limit.

**Proof of the bound.** Regard L as a Laurent series with positive
coefficients, L=L_M+H. If μⱼ=CT Lʲ, coefficient counting gives
μ₄−μ₄,M≤4 CT(HL³). Since μ₂ increases,
c−c_M≤(μ₄−μ₄,M)/24. Also aₖ≤Cpᵏ, with p=√q. Majorize L by C times

```math
P(z)=\sum_{j\in\mathbb Z}p^{|j|}z^{j}
=\frac{1-q}{(1-pz)(1-p/z)}.
```

For k≥0, its cubed coefficient satisfies

```math
[z^{k}]P^{3}=(1-q)^{3}p^{k}
\sum_{r\ge0}\binom{r+k+2}{2}\binom{r+2}{2}q^{r}
\le p^{k}(k+2)^{2}\frac{1+4q+q^{2}}{(1-q)^{2}}
\le 2p^{k}(k+2)^{2}.
```

The first inequality follows from bounding the product of binomials by
(k+2)²(r+1)²(r+2)²/4 and summing its elementary generating function;
the last uses q≤1/20. Sum both signs of k≥N and relax the odd restriction.
The elementary tail ∑_{k≥N}(k+2)²qᵏ gives the proposition.

The exponent improves because a large frequency requires equally large
opposite frequency in a term contributing to CT. This bound does not establish
optimal computational complexity or an all-order monotonicity theorem.

## 6. What is certified

`intervals.py` uses integer endpoints divided by 10¹⁸⁰, with outward rounding
for every operation. It encloses π using Machin's formula and alternating
arctangent remainders, then q using a positive exponential series with a
geometric bound on its tail. Stored digits of π are not an input.

`verify_fourier.py` uses two remainder arguments:

* A uniform wave bound after frequency N−2 is
  ε=8pᴺ/[N(1−qᴺ)(1−q)]. If S₀=2∑ retained aₖ, moment j has error at most
  j(S₀+ε)^{j−1}ε. This verifies logarithmic coefficients of degrees 2,4,6
  using 121 positive harmonics.
* The proposition gives a one-sided quartic enclosure using 61 positive
  harmonics, through frequency 121, with first omitted frequency 123.

Both computations propagate their moment intervals through exact log-series
recurrences. Their reported widths are below 10⁻¹⁵⁰. They agree with an
independent Gauss coefficient recurrence at t=1/2. On |s|≤1/4 its nth raw
term has modulus at most 2^{−n}; thus its degree-j tail after N terms is
at most 4ʲ2^{1−N}, by Cauchy's formula. This supplies a genuine remainder,
not just a comparison at two precisions.

Tests also compare the modular primitive in the reduction audit with the
Gauss value, replay the exact rational witness, and show that an independent-
phase model misses the resonance. Optional mpmath path checks are floating
cross-checks, not certified quadrature. None of these checks is a formal
proof-assistant verification.

## 7. Period normalization and limits

On 0≤x≤r, 1≤T≤α. Replace log(T)ʲ by j ordinary integrals of dt/t over
1≤t≤T. Thus every m₂ⱼ is an absolutely convergent algebraic integral on an
algebraic semialgebraic domain, a Kontsevich–Zagier period [5]. In particular

```math
48K^{2}c=m_0m_4-3m_2^{2}
```

proves K²c is a period. More generally Kᵏℓₖ is a period: the coefficient
of s^{2k} in log(1+∑_{j≥1}m₂ⱼs^{2j}/((2j)!m₀)) has at most k denominator
factors m₀; multiplication by m₀ᵏ leaves a polynomial in periods.
This argument does **not** prove c itself is a period, since arbitrary
inverses are not justified by the ring property, and does not prove
irrationality or independence of these numbers.

The theta, elliptic, Bessel, and period frameworks are established literature.
The application here supplies explicit coordinates, a decomposition, remainder
proofs, and certificates. Priority has not been established. The target is
still not evaluated in a finite familiar-constant basis. The Fourier amplitudes
already contain π, so this is not a standalone formula for computing π.

## Sources

1. NIST DLMF [14.12](https://dlmf.nist.gov/14.12), Legendre integral
   representations, and [15.8](https://dlmf.nist.gov/15.8), hypergeometric
   transformations. Used for the integral representation, not a novelty claim.
2. NIST DLMF [19.8](https://dlmf.nist.gov/19.8), elliptic transformations;
   [22.2](https://dlmf.nist.gov/22.2), Jacobi/theta relations;
   [23.6](https://dlmf.nist.gov/23.6), Weierstrass/Jacobi relations;
   [23.7](https://dlmf.nist.gov/23.7), quarter periods.
3. NIST DLMF [20.5](https://dlmf.nist.gov/20.5), theta products. The logarithmic
   Fourier expansion here is derived from those products with branches fixed.
4. NIST DLMF [10.35.1](https://dlmf.nist.gov/10.35.E1), modified Bessel
   generating function.
5. M. Kontsevich and D. Zagier, *Periods*, Mathematics Unlimited—2001 and Beyond
   (2001), 771–808. [Author PDF](https://www.ihes.fr/~maxim/TEXTS/Periods.pdf).

More specialized sources, including torsion-point elliptic polylogarithms and
boundary constants, are discussed in [BOUNDED_REDUCTION.md](BOUNDED_REDUCTION.md).
