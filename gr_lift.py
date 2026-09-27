#!/usr/bin/env python3
"""Evaluate a coefficientwise lift directly on GR(4, m), without the derivative criterion.

GR(4, m) = Z_4[xi]/(h), where h is the modulus of the table read with 0/1 coefficients in Z_4.
The lift of the reduced representative f = sum a_k x^k replaces each coefficient a_k, whose
bits are its coordinates in the basis 1, alpha, ..., alpha^(m-1), by the element of GR(4, m)
with the same 0/1 coordinates in 1, xi, ..., xi^(m-1). The script evaluates this lift at all
4^m ring elements and counts distinct values. A permutation of GR(4, m) refutes Conjecture 1 of
arXiv:2608.30808 directly, since k = 2 is the case that every k > 1 reduces to.

    python3 gr_lift.py          # hand example, apn_m5, apn_m7 and the x^3 control; seconds
    python3 gr_lift.py --all    # adds apn_m9 (262,144 ring elements); about 10 minutes

Exits non-zero if a lift of a counterexample is not a permutation, or if the control is.
"""
import itertools
import json
import os
import sys

import check


def gr_mul(a, b, h, m):
    """a, b: m coordinates mod 4; h: the low m coefficients of the monic modulus."""
    prod = [0] * (2 * m - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    prod[i + j] = (prod[i + j] + ai * bj) & 3
    for k in range(2 * m - 2, m - 1, -1):  # xi^m = -sum h_i xi^i
        c = prod[k]
        if c:
            prod[k] = 0
            for i in range(m):
                if h[i]:
                    prod[k - m + i] = (prod[k - m + i] - c * h[i]) & 3
    return prod[:m]


def lift_image_size(m, mod, T):
    coef = check.reduced_representative(m, mod, T)
    h = [(mod >> i) & 1 for i in range(m)]
    lift = [[(c >> i) & 1 for i in range(m)] for c in coef]
    deg = max(k for k in range(1 << m) if coef[k])
    seen = set()
    for z in itertools.product(range(4), repeat=m):
        z = list(z)
        acc = [0] * m
        for k in range(deg, -1, -1):
            acc = gr_mul(acc, z, h, m)
            acc = [(u + v) & 3 for u, v in zip(acc, lift[k])]
        seen.add(tuple(acc))
    return len(seen)


def main():
    failures = 0

    def report(label, m, mod, T, want_permutation):
        nonlocal failures
        n = lift_image_size(m, mod, T)
        ok = (n == 4 ** m) == want_permutation
        failures += not ok
        print('  %s  %s: %d distinct values on the %d elements of GR(4, %d)'
              % ('ok  ' if ok else 'FAIL', label, n, 4 ** m, m), flush=True)

    m, mod = 5, 37
    L = lambda x: x ^ check.power(x, 8, m, mod) ^ check.power(x, 16, m, mod)
    report('hand example (x + x^8 + x^16)^3, lift permutes', m, mod,
           [check.power(L(x), 3, m, mod) for x in range(32)], True)
    names = ['apn_m5.json', 'apn_m7.json'] + (['apn_m9.json'] if '--all' in sys.argv else [])
    d = check.tables_dir()
    for name in names:
        doc = json.load(open(os.path.join(d, name)))
        report(name + ', lift permutes', doc['m'], doc['modulus_bits'], doc['table'], True)
    # control: x^3 is an APN permutation of GF(8) with a critical point, so its lift must not permute
    report('control x^3 over GF(8), lift does not permute', 3, 0b1011,
           [check.power(x, 3, 3, 0b1011) for x in range(8)], False)
    print('all lifts behave as stated' if not failures else '%d lift check(s) FAILED' % failures)
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
