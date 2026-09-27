# The Reduced Critical-Point Conjecture is false

Bartoli and Stănică, *Reduced polynomial lifts of APN permutations over Galois rings and effective non-APN bounds* ([arXiv:2608.30808v1](https://arxiv.org/abs/2608.30808), 31 August 2026), state in their LaTeX source (lines 161–164):

```latex
\begin{conj}[Reduced Critical-Point Conjecture]
\label{conj:critical}
For every $q=2^m$, the reduced representative $f\in\F_q[x]$ of every APN permutation of $\F_q$ has a critical point in $\F_q$.
\end{conj}
```

The reduced representative is the unique polynomial of degree less than `q` that induces the function.
A critical point is a zero of the formal derivative `f'`.
By the paper's Theorem 2.2, this statement is equivalent to its Conjecture 1 (lines 113–116):

```latex
\begin{conj}[Reduced APN Lifting Conjecture]
\label{conj:main}
Let $q=2^m$, and let $f\in\F_q[x]$ be the reduced representative, of degree less than $q$, of an APN permutation of $\F_q$. Then no coefficientwise lift of $f$ induces a permutation of $\GR(2^k,m)$ for any $k>1$.
\end{conj}
```

Conjecture 1 is the normalized form of Conjecture 2 of Rønjom and Sandrib, *On the differential uniformity of polynomials over Galois rings* ([Cryptography and Communications, 2026](https://doi.org/10.1007/s12095-026-00895-x); [IACR ePrint 2026/1393](https://eprint.iacr.org/2026/1393)).

**This repository refutes both statements, already at m = 5.**

## A counterexample you can check by hand

Over `F_32`, take

```
F(x) = (x + x^8 + x^16)^3.
```

1. `L(x) = x + x^8 + x^16` is `F_2`-linear.
   It is a bijection of `F_32`, because its associated polynomial `1 + t^3 + t^4` is coprime to `t^5 + 1` over `F_2`.
   (`t^5 + 1 = (t + 1)(t^4 + t^3 + t^2 + t + 1)`, and `t^4 + t^3 + 1` is irreducible and is neither factor.)
2. `x^3` is an APN permutation of `F_32` (a Gold exponent, and `gcd(3, 31) = 1`).
   A linear bijection on the input keeps both properties, so `F` is an APN permutation.
3. On `F_32`, `x^32 = x`, so `(x + x^8 + x^16)^2 = x^2 + x^16 + x`. Then

   ```
   F(x) = (x + x^2 + x^16)(x + x^8 + x^16)
        = x^2 + x^9 + x^17 + x^3 + x^10 + x^18 + x^17 + x^24 + x^32.
   ```

   The two `x^17` terms cancel and `x^32 = x`, so the reduced representative is

   ```
   f(x) = x^24 + x^18 + x^10 + x^9 + x^3 + x^2 + x.
   ```

   The derivation uses no modulus, so the same `f` holds for every choice of modulus for `F_32`.
4. In characteristic 2 only the odd exponents survive differentiation:

   ```
   f'(x) = x^8 + x^2 + 1 = (x^4 + x + 1)^2.
   ```

5. `x^4 + x + 1` is irreducible over `F_2`, so its roots lie in `F_16` and not in `F_2`.
   Since `F_16 ∩ F_32 = F_2`, `f'` has no zero in `F_32`.
   So `f` has no critical point, and Conjecture 2 fails.
6. By the paper's Theorem 2.1 (the Galois-ring permutation criterion), every coefficientwise lift of `f` permutes `GR(2^k, 5)` for every `k > 1`.
   So Conjecture 1 fails too.
   `gr_lift.py` confirms this without the criterion: one lift takes 1024 distinct values on the 1024 elements of `GR(4, 5)`.

## Why it happens

Section 7 of the paper covers power maps: for `f(x) = a·x^e` with `e` odd, `f'(0) = 0`.
That argument does not survive an `F_2`-linear change of variable on the input.
The example above is `x^3` composed with the linear bijection `L`.
Having a critical point is a property of one representative, not of an equivalence class.

Two kinds of relabelling always leave a critical point:

- An `F_2`-linear map on the output side alone.
  `B(x^d)` is a sum of terms `x^(d·2^i)`, and no exponent reduces to 1 unless `x^d` is itself linear.
  So the reduced representative has no `x` term, and `f'(0) = 0`.
- An `F_q`-affine map `x ↦ ax + b` on either side.
  On the output side it leaves the zeros of `f'` unchanged.
  On the input side it moves the critical point of a power map from `0` to `b/a`.

Relabelling the input, or both sides, often removes every critical point.
Random two-sided `F_2`-affine relabellings of every APN power map had no critical point about 37% of the time:
370,430 of 1,000,000 at m = 5, 817,835 of 2,200,000 at m = 7, and 242,811 of 650,000 at m = 9.
Input-side relabellings alone did so 1,540,572 times in 3,850,000 samples (40%).
Output-side relabellings and `F_q`-affine relabellings did so 0 times in 3,850,000 samples each.
These samples come from the audit of this result; `verify.py` does not rerun them.

## The smallest field

At m = 3 the conjecture holds.
All 40,320 permutations of `F_8` include 10,752 APN permutations, and each has a critical point.
`check.py` verifies this exhaustively.
`F_4` and `F_16` have no APN permutations.
So m = 5 is the smallest field with a counterexample, setting aside `F_2`, where every map is affine.

## Three more examples

Each `apn_m*.json` file gives a modulus and a value table, `table[i] = F(i)`.
The integer `i` stands for the field element `sum of b_j·α^j`, where `b_j` is bit `j` of `i` and `α` is a root of the modulus.

| file | field and modulus | reduced degree | critical points | `F_2`-affine equivalent to |
| --- | --- | --- | --- | --- |
| `apn_m5.json` | `F_32`, `x^5 + x^2 + 1` | 28 | 0 | `x^11` (Kasami class) |
| `apn_m7.json` | `F_128`, `x^7 + x + 1` | 112 | 0 | `x^11` (Welch) |
| `apn_m9.json` | `F_512`, `x^9 + x^4 + 1` | 384 | 0 | `x^5` (Gold) |

Each is an APN permutation.
`check.py` also re-evaluates, at every point, a certificate of the form `F(x) = B((A x + c)^d) + t` with `A` and `B` invertible `F_2`-linear maps.
The lifts to `GR(4, 5)`, `GR(4, 7)` and `GR(4, 9)` take 1024, 16,384 and 262,144 distinct values: they are permutations.

## Verify it yourself

```sh
python3 verify.py          # a few seconds
python3 verify.py --all    # adds the GR(4, 9) lift; about 10 minutes in pure Python
```

Python 3.8 or later, standard library only.
Each command exits non-zero if any check fails.
CI runs `python3 verify.py --all` on every push.

- `check.py` builds `F_(2^m)` arithmetic from scratch.
  For each table and for the hand example, it tests that the modulus is irreducible, that `F` is a permutation with differential uniformity 2, and that `f'` has no zero.
  It computes `f` by Lagrange interpolation and re-evaluates it at every point, so the coefficients certify themselves.
  It also runs the m = 3 exhaustion.
- `gr_lift.py` evaluates a coefficientwise lift directly on all of `GR(4, m)` and counts distinct values.
  As a control, it checks that the lift of `x^3` over `F_8`, which has a critical point, does not permute `GR(4, 3)`.
- `verify.py` runs both.

## What this is, and what it is not

- **One result, not three.** The counterexamples are elementary.
  Each table is a relabelled textbook power map, so none is a new APN function.
  The three tables and the hand example show one mechanism.
- **Not claimed as new:** that critical points are not invariant under `F_2`-affine equivalence.
- **Claimed:** Conjecture 2 of arXiv:2608.30808, as stated there, is false, with explicit witnesses at m = 5, 7 and 9.
  So are its Conjecture 1 and the normalized form of the Rønjom–Sandrib conjecture.
- **Novelty:** a search on 2026-09-27 of arXiv, IACR ePrint, Semantic Scholar, OpenAlex and the web found no earlier counterexample.
  Work that is recent, unindexed or private cannot be excluded.
- **Open:** we have not studied corrected forms of the conjecture.

## Provenance

The three tables were found with Tangle's automated research system, by a search over random `F_2`-affine relabellings of APN power maps.
They were then checked independently, by separately written programs that share no code and each exhaust the field.
`check.py` is one of those programs.
The hand example came from the audit of the mechanism.

## Author

Drew Stone — Tangle — drew@tangle.tools

## License

MIT.
