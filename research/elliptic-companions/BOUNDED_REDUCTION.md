# A bounded reduction attempt: result and limits

**Result:** no finite evaluation of the special value c was obtained. The
attempt fixes a canonical torsion alphabet, shows exactly why a naive
quarter-turn does not close on one logarithmic differential, reduces the
independent-phase contribution to ordinary Eisenstein primitives, and excludes
a stated 28-function functional ansatz by exact rational arithmetic.

This is a stopping point for this investigation, not an impossibility theorem
for the number c. In particular, no complete Fay-relation calculation or
minimal-basis theorem for all torsion-point iterated integrals is claimed.

## 1. A concrete torsion alphabet, with the path retained

Use the curve and conventions of [NOTE.md](NOTE.md). In normalized torus
coordinates u=z/(2ω), γ runs from u₁=(−1+i)/4 to u₀=i/4. Put v=u−u₁,
so the new path is the real interval [0,1/4]. Define

```math
f_1(t)=\pi\frac{\theta_1'(\pi t,e^{-\pi})}
{\theta_1(\pi t,e^{-\pi})}.
```

The prime differentiates the theta function's first argument. The divisor of g
and the theta products give the exact differential

```math
\eta=d\log g=\sum_{j=1}^{4}\epsilon_j f_1(v-b_j)\,dv,
\qquad (\epsilon_1,\epsilon_2,\epsilon_3,\epsilon_4)=(1,1,-1,-1),
```
```math
(b_1,b_2,b_3,b_4)=
\left(\frac{3-i}{4},\frac{3+i}{4},
\frac{1+i}{4},\frac{1-i}{4}\right).
```

These are specified lifts, not silently reduced modulo periods. Their signed
sum is 1, and their signed imaginary parts sum to zero. Consequently the
ratio of theta products is elliptic with the same divisor as g and differs
by a constant; its logarithmic derivative is exactly η. There is no missing
holomorphic differential. The completed, doubly periodic first Kronecker
kernel gives the same signed sum because its nonholomorphic corrections cancel.
None of the four poles lies on the real integration segment.

For the iterated-integral convention with the leftmost differential outermost,
let Jⱼ be the sum of words (dv, f₁(v−bₖ₁)dv, …, f₁(v−bₖⱼ)dv), integrated
from 0 to 1/4 and weighted by εₖ₁⋯εₖⱼ. Then Iⱼ=−2ωJⱼ and J₀=1/4, so

```math
R=16J_2,\qquad c=64J_4-128J_2^{2}.
```

This is an explicit expression in the published torsion-point elliptic
iterated-integral framework [1,2]: 16 words for J₂ and 256 for J₄ before
relations. These counts are expression sizes, **not dimensions or minimal
basis counts**. The conversion fixes the alphabet and boundary, but does
not itself reduce the number of independent unknown special values.

## 2. Why quarter-turn symmetry does not evaluate the integral by itself

The square curve's rotation is ρ(X,Y)=(−X,iY). For g=−iY/(X+1), direct
substitution using Y²=X³−X gives

```math
g(\rho(P))g(P)=iX(P).
```

Thus, setting ξ=dlog X,

```math
\rho^{\ast}\eta=\xi-\eta,\qquad
\rho^{\ast}\xi=\xi.
```

The extra differential is not a constant multiple of η. The residue vectors
at (1,0), (0,0), (−1,0), O are respectively

| Differential | (1,0) | (0,0) | (−1,0) | O |
|---|---:|---:|---:|---:|
| η | 1 | 1 | −1 | −1 |
| ξ | 0 | 2 | 0 | −2 |

They have rank two. These pole residues establish the failure of one-letter
closure without numerical fitting. The rotation also sends γ to iγ, changing
both endpoints and the lift. Taking logarithms of the product identity
requires continuing the branch along that path; replacing it by a real
logarithmic identity with no additive constant is unjustified.

Translations by (1,0) and (−1,0) give g↦1/g and g↦−g, respectively. The
first recovers sign symmetry of logarithmic moments. The second preserves
dlog g but can change the logarithm by a branch constant. These transformations
do not by themselves remove J₄. `audit_geometry.py` checks all three actions
symbolically. This does not rule out a more elaborate relation in a larger
alphabet.

## 3. Separating the ordinary part of the resonance formula

For the following functional test let q vary, D=q d/dq, and let D^{−1}
mean the primitive with zero constant coefficient. Define

```math
A_r(q)=\sum_{n\ge1}
\bigl(\sigma_r(n)-\sigma_r(n/2)\bigr)q^{n},
\qquad V(q)=16A_3(q),\qquad R(q)=D^{-2}V(q).
```

Here σᵣ(n/2)=0 for odd n. At q=e^{−π}, R(q) agrees with R in the note.
The modular quartic extension is

```math
C(q)=-D^{-2}\bigl[(D^{-1}V(q))^{2}\bigr],\qquad C(e^{-\pi})=c.
```

To justify the last equality, use t=λ(τ), q=e^{πiτ}, F₀=F(1/2,1/2;1;t)
and D t=t(1−t)F₀². Substituting
h(s,q)=(1−t)ˢ F(1/2+s,1/2+s;1;t)/F₀ into Gauss's equation gives
D²h=s²tF₀⁴h, h(s,0)=1. The classical theta identity
V=θ₂⁴θ₃⁴=16A₃ identifies the kernel. Equating the quadratic and quartic
coefficients of log h proves the displayed formulas. At τ=i, t=1/2.
The theta identity also follows by combining the duplication formulas with
E₄(τ)=(θ₂⁸+θ₃⁸+θ₄⁸)/2: explicitly
V=(E₄(τ/2)−E₄(τ))/15, with the usual E₄ Fourier convention.

The fourth-power term in the independent-phase model reduces exactly:

```math
\sum_{k>0,\ k\ \mathrm{odd}}a_k(q)^{4}
=\frac{128}{3}D^{-4}\bigl(A_7(q)-A_5(q)\bigr).
```

Proof: expand q^{2k}/(1−qᵏ)⁴ as
∑_{m≥1}(m³−m)q^{km}/6 and reindex n=km. This is an analytic identity for
|q|<1, not an inference from a finite coefficient check. The script separately
checks the first 60 coefficients.

For explicit Eisenstein normalization, if
𝔈₆(q)=1−504∑σ₅(n)qⁿ and 𝔈₈(q)=1+480∑σ₇(n)qⁿ, the right side is

```math
\frac{4}{45}D^{-4}\bigl(\mathcal E_8(q)-\mathcal E_8(q^{2})\bigr)
+\frac{16}{189}D^{-4}\bigl(\mathcal E_6(q)-\mathcal E_6(q^{2})\bigr).
```

This removes the four-amplitude sum from that part of the formula, but adds
no evaluation of the remaining interaction. A single ordinary primitive is
not automatically a known constant at this endpoint.

## 4. An exact, explicitly limited exclusion certificate

Test membership of C(q) in the rational linear span of these **28 functions**:

* 1;
* D^{−j}Aᵣ(qᵈ), for r∈{3,5,7}, d∈{1,2}, j∈{1,2,3,4};
* R(q)², R(q)R(q²), R(q²)².

This family contains the preceding ordinary fourth-power contribution, ordinary
primitives at two scales and three weights, and the natural quadratic products.
It does not include arbitrary algebraic coefficient functions or every
possible modular form or iterated integral.

`verify_bounded_reduction.py` computes exact rational coefficients through q⁶⁰.
The candidate matrix has rank **28**; adjoining C gives rank **29**. It also
constructs a rational linear functional on coefficients, supported on odd
powers q¹,q³,…,q²⁷, that is zero on each proposed basis function and exactly
one on C. The complete rational weights are in
[the certificate](certificates/bounded_reduction.json).

This is a proof of nonmembership in the stated **functional** span: any
identity of full power series would have to satisfy the witnessed finite
coefficient identity. It is stronger than a failed floating-point search
within that ansatz, but it does **not** exclude:

* a different or larger function family;
* nonlinear identities or different coefficient fields;
* a relation that holds only at q=e^{−π};
* an eventual reduction using full torsion-point Fay relations and boundary data.

The regression suite reconstructs every candidate vector and checks the
witness directly, independently of its row-reduction construction. A tampered
target coefficient is a negative control. Another test checks C(e^{−π})
against the Gauss-series value with a rigorous tail: Vₙ≤32n³ implies
|(C)ₙ|≤35n³. After N terms use
35N³qᴺ/[1−q(1+1/N)³]. No experimental failure is presented as numerical
independence of constants.

## 5. Literature assessment and stop decision

Broedel et al. [1] provide the torsion-point elliptic-polylogarithm to
Eisenstein-integral framework. Their section 8.2 is a worked conversion,
not a universal theorem that every such special value reduces to ordinary
zeta or gamma values. Their discussion also explains why a symbol alone
loses constants.

Broedel–Matthes–Richter–Schlotterer [2] treat twisted elliptic multiple zeta
values and Fay relations. The translated quarter-period alphabet here fits
that general setting, but we have not implemented or exhausted all those
relations for the 256-word expression.

Duhr–Tancredi [3], sections 3–5, show explicitly how changes of variables
and modular transformations require boundary/regularization constants.
Our finite path avoids singular endpoints, but subsequent conversion to
cusp-based Eisenstein integrals must still keep its boundary data.

Rogers–Wan–Zucker [4] prove closed evaluations for certain moments of complete
elliptic integrals. Those are useful examples of successful reductions, not
an evaluation of this logarithmic theta-ratio moment.

The framework is established. The contribution recorded here is an explicit
application, exact symmetry and finite-ansatz audits, and reproducible
certification. This sweep establishes no priority claim. We stop with c
unreduced; further progress should eliminate a required special value or
prove a new relation, rather than merely rename an infinite sum.

## Sources

1. J. Broedel, C. Duhr, F. Dulat, B. Penante, L. Tancredi,
   *Elliptic symbol calculus: from elliptic polylogarithms to iterated
   integrals of Eisenstein series*, JHEP 08 (2018) 014, §§3,8–9.
   [arXiv:1803.10256](https://arxiv.org/abs/1803.10256).
2. J. Broedel, N. Matthes, G. Richter, O. Schlotterer,
   *Twisted elliptic multiple zeta values and non-planar one-loop open-string
   amplitudes*, J. Phys. A 51 (2018) 285401, §3 and appendix B.
   [arXiv:1704.03449](https://arxiv.org/abs/1704.03449).
3. C. Duhr, L. Tancredi, *Algorithms and tools for iterated Eisenstein
   integrals*, JHEP 02 (2020) 105, §§3–5.
   [arXiv:1912.00077](https://arxiv.org/abs/1912.00077).
4. M. Rogers, J. G. Wan, I. J. Zucker, *Moments of elliptic integrals and
   critical L-values* (2013 preprint).
   [arXiv:1303.2259](https://arxiv.org/abs/1303.2259).
