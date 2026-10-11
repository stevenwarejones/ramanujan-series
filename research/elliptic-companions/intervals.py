"""Outward integer intervals and an independent Gauss-series reference.

All endpoints are integers divided by 10**180. Only the standard library
is needed. No stored decimal expansion of any special constant is used.
"""
from dataclasses import dataclass
from math import isqrt

PLACES = 180
SCALE = 10**PLACES


def ceil_div(a, b):
    assert b > 0
    return -((-a) // b)


@dataclass(frozen=True)
class I:
    lo: int
    hi: int

    def __post_init__(self):
        assert self.lo <= self.hi

    @staticmethod
    def exact(n):
        return I(n * SCALE, n * SCALE)

    def __add__(self, other):
        other = as_i(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-as_i(other))

    def __rsub__(self, other):
        return as_i(other) - self

    def __mul__(self, other):
        other = as_i(other)
        products = [x * y for x in (self.lo, self.hi)
                    for y in (other.lo, other.hi)]
        return I(min(products) // SCALE, ceil_div(max(products), SCALE))

    __rmul__ = __mul__

    def reciprocal(self):
        assert self.lo > 0, 'Only positive denominators are needed here.'
        return I(SCALE * SCALE // self.hi,
                 ceil_div(SCALE * SCALE, self.lo))

    def __truediv__(self, other):
        return self * as_i(other).reciprocal()

    def __rtruediv__(self, other):
        return as_i(other) / self

    def __pow__(self, exponent):
        assert isinstance(exponent, int) and exponent >= 0
        result = I.exact(1)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent >>= 1
        return result

    def sqrt(self):
        assert self.lo >= 0
        lower = isqrt(self.lo * SCALE)
        upper = isqrt(self.hi * SCALE)
        if upper * upper < self.hi * SCALE:
            upper += 1
        return I(lower, upper)



def as_i(value):
    return value if isinstance(value, I) else I.exact(value)


def convolution(a, b):
    assert len(a) == len(b)
    return [sum((a[k] * b[j-k] for k in range(j+1)), I.exact(0))
            for j in range(len(a))]



def decimal_string(n):
    sign = '-' if n < 0 else ''
    n = abs(n)
    return f'{sign}{n // SCALE}.{n % SCALE:0{PLACES}d}'


def pi_interval(terms=150):
    def atan_inverse(m):
        x = I.exact(1)/m
        power = x
        total = I.exact(0)
        for n in range(terms):
            total += ((-1)**n)*power/(2*n+1)
            power *= x*x
        tail = power/(2*terms+1)
        return I(total.lo-tail.hi, total.hi+tail.hi)
    return 16*atan_inverse(5)-4*atan_inverse(239)


def exp_positive(value, terms=280):
    assert 0 <= value.lo <= value.hi < 4*SCALE
    term = I.exact(1)
    total = I.exact(0)
    for n in range(terms):
        total += term
        term = term*value/(n+1)
    tail = term/(1-value/(terms+1))
    return I(total.lo, total.hi+tail.hi)


def log_jet(raw):
    out = [I.exact(0)]*len(raw)
    for n in range(1,len(raw)):
        out[n] = (n*raw[n]-sum(
            (k*out[k]*raw[n-k] for k in range(1,n)),I.exact(0)))/(n*raw[0])
    return out


def gauss_log_jet(terms=650, order=6):
    """Log F(1/2+s,1/2+s;1;1/2); two-sided Cauchy tail on |s|<=1/4."""
    zero = I.exact(0)
    current = [I.exact(1)]+[zero]*order
    total = [zero]*(order+1)
    for n in range(terms):
        total = [x+y for x,y in zip(total,current)]
        numerator = [I.exact((2*n+1)**2)/4,I.exact(2*n+1),I.exact(1)]+[zero]*(order-2)
        current = [x/(2*(n+1)**2) for x in convolution(current,numerator)]
    for j in range(order+1):
        tail = (I.exact(4**j)*2/2**terms).hi
        total[j] = I(total[j].lo-tail,total[j].hi+tail)
    return log_jet(total)
