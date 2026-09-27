#!/usr/bin/env python3
"""Check the counterexamples to the Reduced Critical-Point Conjecture (arXiv:2608.30808, Conjecture 2).

Python 3.8+, standard library only, a few seconds. Exits non-zero if any check fails.

For each value table apn_m5.json, apn_m7.json, apn_m9.json beside this file, and for the
hand-checkable example F(x) = (x + x^8 + x^16)^3 over GF(32), it checks:
  - the modulus is irreducible (Rabin's test), so the table lives on a field;
  - F is a permutation and its differential uniformity is 2 (APN), over every pair (a, x);
  - the reduced representative f (degree < q), computed by Lagrange interpolation, reproduces
    F at every point, so its coefficients certify themselves;
  - its formal derivative f' has no zero in GF(2^m).
It also re-evaluates, at every point, a certificate that each table is F_2-affine equivalent
to a power map, and checks the smallest field exhaustively: at m = 3 every APN permutation
of GF(8) has a critical point.

Field elements are integers: bit j of i is the coefficient of alpha^j, alpha a root of the modulus.
"""
import itertools
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def mul(a, b, m, mod):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> m:
            a ^= mod
    return r


def power(a, e, m, mod):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a, m, mod)
        a = mul(a, a, m, mod)
        e >>= 1
    return r


def pmod(a, b):
    db = b.bit_length() - 1
    while a and a.bit_length() - 1 >= db:
        a ^= b << (a.bit_length() - 1 - db)
    return a


def pgcd(a, b):
    while b:
        a, b = b, pmod(a, b)
    return a


def pmulmod(a, b, mod):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a.bit_length() > mod.bit_length() - 1:
            a = pmod(a, mod)
    return pmod(r, mod)


def rabin_irreducible(f):
    """f of degree m is irreducible iff x^(2^m) = x mod f and gcd(x^(2^(m/p)) - x, f) = 1
    for each prime p dividing m."""
    m = f.bit_length() - 1
    x = 2

    def frob(k):
        y = x
        for _ in range(k):
            y = pmulmod(y, y, f)
        return y
    if frob(m) != pmod(x, f):
        return False
    primes = [p for p in range(2, m + 1) if m % p == 0 and all(p % d for d in range(2, p))]
    return all(pgcd(f, frob(m // p) ^ x) == 1 for p in primes)


def reduced_representative(m, mod, T):
    """Coefficients a_0..a_(q-1) of the unique polynomial of degree < q with f(x) = T[x].
    a_0 = F(0); a_j = sum_{c != 0} F(c) c^(q-1-j) for 1 <= j <= q-2; a_(q-1) = sum_c F(c)."""
    q = 1 << m
    coef = [0] * q
    coef[0] = T[0]
    for v in T:
        coef[q - 1] ^= v
    for c in range(1, q):
        if T[c] == 0:
            continue
        p = 1
        for j in range(q - 2, 0, -1):
            p = mul(p, c, m, mod)  # c^(q-1-j)
            coef[j] ^= mul(T[c], p, m, mod)
    return coef


def horner(cs, x, m, mod):
    r = 0
    for c in reversed(cs):
        r = mul(r, x, m, mod) ^ c
    return r


def differential_uniformity(T):
    q = len(T)
    du = 0
    for a in range(1, q):
        counts = [0] * q
        for x in range(q):
            counts[T[x ^ a] ^ T[x]] += 1
        du = max(du, max(counts))
    return du


def analyse(m, mod, T):
    q = 1 << m
    assert len(T) == q and all(0 <= v < q for v in T)
    coef = reduced_representative(m, mod, T)
    deriv = [coef[j + 1] if (j + 1) % 2 == 1 else 0 for j in range(q - 1)]  # f' in characteristic 2
    return {
        'modulus_irreducible': rabin_irreducible(mod),
        'permutation': len(set(T)) == q,
        'differential_uniformity': differential_uniformity(T),
        'reduced_rep_reproduces_table': all(horner(coef, x, m, mod) == T[x] for x in range(q)),
        'reduced_degree': max((j for j in range(q) if coef[j]), default=-1),
        'critical_points': sum(1 for x in range(q) if horner(deriv, x, m, mod) == 0),
    }, coef


def linear(images, v):
    r = 0
    i = 0
    while v:
        if v & 1:
            r ^= images[i]
        v >>= 1
        i += 1
    return r


def is_equivalent_to_power(m, mod, T, d, A, c2, B, t1):
    """T[x] == B((A x + c2)^d) + t1 at every x, with A and B invertible F_2-linear maps
    (given by the images of the basis vectors 1, alpha, ..., alpha^(m-1))."""
    q = 1 << m
    return (len({linear(A, v) for v in range(q)}) == q
            and len({linear(B, v) for v in range(q)}) == q
            and all(T[x] == linear(B, power(linear(A, x) ^ c2, d, m, mod)) ^ t1 for x in range(q)))


# name: (expected reduced degree, power d, A, c2, B, t1)
TABLES = {
    'apn_m5.json': (28, 11, [1, 2, 20, 4, 24], 0, [6, 11, 2, 18, 26], 23),
    'apn_m7.json': (112, 11, [25, 24, 81, 13, 70, 105, 76], 1, [70, 26, 119, 126, 2, 88, 54], 110),
    'apn_m9.json': (384, 5, [66, 177, 209, 17, 327, 345, 490, 287, 316], 1,
                    [198, 98, 413, 300, 321, 470, 167, 437, 355], 365),
}

FAILURES = []


def expect(label, ok):
    print(('  ok    ' if ok else '  FAIL  ') + label)
    if not ok:
        FAILURES.append(label)


def expect_counterexample(res, degree):
    expect('modulus irreducible', res['modulus_irreducible'])
    expect('permutation', res['permutation'])
    expect('APN (differential uniformity %d)' % res['differential_uniformity'], res['differential_uniformity'] == 2)
    expect('reduced representative reproduces the table at every point', res['reduced_rep_reproduces_table'])
    expect('reduced degree %d (expected %d)' % (res['reduced_degree'], degree), res['reduced_degree'] == degree)
    expect("f' has %d zeros in the field (a counterexample needs 0)" % res['critical_points'], res['critical_points'] == 0)


def tables_dir():
    for d in (HERE, os.path.join(HERE, 'tables')):
        if all(os.path.exists(os.path.join(d, n)) for n in TABLES):
            return d
    sys.exit('apn_m5.json, apn_m7.json and apn_m9.json must sit beside check.py')


def check_tables():
    d = tables_dir()
    for name, (degree, pw, A, c2, B, t1) in TABLES.items():
        doc = json.load(open(os.path.join(d, name)))
        m, mod, T = doc['m'], doc['modulus_bits'], doc['table']
        print('%s: GF(2^%d), modulus %s' % (name, m, doc['modulus']))
        res, _ = analyse(m, mod, T)
        expect_counterexample(res, degree)
        expect('F_2-affine equivalent to x^%d (certificate checked at every point)' % pw,
               is_equivalent_to_power(m, mod, T, pw, A, c2, B, t1))


def check_hand_example():
    """F(x) = (x + x^8 + x^16)^3 over GF(32), under every irreducible quintic modulus."""
    print('hand example F(x) = (x + x^8 + x^16)^3 over GF(32), under each irreducible quintic modulus')
    m = 5
    moduli = [f for f in range(32, 64) if rabin_irreducible(f)]
    for mod in moduli:
        L = lambda x: x ^ power(x, 8, m, mod) ^ power(x, 16, m, mod)
        T = [power(L(x), 3, m, mod) for x in range(32)]
        res, coef = analyse(m, mod, T)
        ok = (len({L(x) for x in range(32)}) == 32
              and res['permutation'] and res['differential_uniformity'] == 2
              and res['reduced_rep_reproduces_table'] and res['critical_points'] == 0
              and [j for j in range(32) if coef[j]] == [1, 2, 3, 9, 10, 18, 24]
              and all(coef[j] == 1 for j in range(32) if coef[j]))
        expect('modulus bits %d: L is a bijection, F is an APN permutation, '
               'f = x^24 + x^18 + x^10 + x^9 + x^3 + x^2 + x, f\' has no zero' % mod, ok)
    expect('%d irreducible quintic moduli (expected 6)' % len(moduli), len(moduli) == 6)


def check_m3():
    """Every permutation of GF(8): count the APN ones and check each has a critical point."""
    m, mod = 3, 0b1011
    apn = with_critical_point = 0
    for T in itertools.permutations(range(8)):
        if differential_uniformity(T) != 2:
            continue
        apn += 1
        coef = reduced_representative(m, mod, T)
        deriv = [coef[j + 1] if (j + 1) % 2 == 1 else 0 for j in range(7)]
        if any(horner(deriv, x, m, mod) == 0 for x in range(8)):
            with_critical_point += 1
    print('GF(8), modulus x^3 + x + 1: all 40320 permutations')
    expect('%d APN permutations (expected 10752)' % apn, apn == 10752)
    expect('%d of them have a critical point (the conjecture holds at m = 3)' % with_critical_point,
           with_critical_point == apn)


def main():
    check_tables()
    check_hand_example()
    check_m3()
    if FAILURES:
        print('\n%d check(s) FAILED' % len(FAILURES))
        return 1
    print('\nall checks passed: each table and the hand example is an APN permutation whose reduced '
          'representative has no critical point, so Conjecture 2 of arXiv:2608.30808 is false')
    return 0


if __name__ == '__main__':
    sys.exit(main())
