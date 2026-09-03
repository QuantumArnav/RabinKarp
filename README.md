# The Karp–Rabin Fingerprint Algorithm — Derivation, Collision Analysis, and Implementation

## 1. The problem

Given a text $T$ of length $n$ and a pattern $P$ of length $m$ (both over an
alphabet of size $\sigma$, e.g. $\sigma = 256$ for bytes), find every start
index $i$ such that $T[i \dots i+m-1] = P$.

The naive algorithm compares $P$ against every window of $T$ character by
character: $O(nm)$ worst case. Karp and Rabin's idea (1987; building on
Rabin's 1981 fingerprinting technique) is to give every length-$m$ window a
cheap numeric **fingerprint**, compare fingerprints in $O(1)$, and only fall
back to a full character comparison when two fingerprints agree.

## 2. Representing a window as a number

Fix a base $b \ge \sigma$ (so every character maps to a distinct "digit")
and a modulus $q$. Treat the $m$ characters of a window as the digits of a
base-$b$ integer, reduced mod $q$:

$$
h(T[i \dots i+m-1]) \;=\; \Big(\sum_{k=0}^{m-1} T[i+k]\cdot b^{\,m-1-k}\Big) \bmod q
$$

Computed directly with Horner's rule, this costs $O(m)$ per window — no
better than brute force if we recompute it from scratch at every $i$. The
whole point of the algorithm is to update it in $O(1)$ as the window slides.

## 3. Deriving the rolling-hash recurrence

Write $h_i$ for the fingerprint of the window starting at $i$:

$$
h_i = T[i]\,b^{m-1} + T[i+1]\,b^{m-2} + \dots + T[i+m-1]\,b^{0} \pmod q
$$

$$
h_{i+1} = T[i+1]\,b^{m-1} + T[i+2]\,b^{m-2} + \dots + T[i+m]\,b^{0} \pmod q
$$

Multiply $h_i$ by $b$:

$$
b\cdot h_i = T[i]\,b^{m} + T[i+1]\,b^{m-1} + \dots + T[i+m-1]\,b^{1}
$$

Compare this to $h_{i+1}$: every term of $b\cdot h_i$ except the first,
$T[i]\,b^m$, is exactly one power of $b$ higher than the matching term of
$h_{i+1}$'s... more directly, subtract the *leading* digit's contribution
from $h_i$ **before** multiplying, then add the new trailing digit:

$$
h_{i+1} = \Big(\big(h_i - T[i]\cdot b^{m-1}\big)\cdot b + T[i+m]\Big) \bmod q
$$

**Check:** expand the right-hand side.
$\big(h_i - T[i]b^{m-1}\big)$ removes $T[i]$'s term, leaving
$T[i+1]b^{m-2}+\dots+T[i+m-1]b^0$ (as a value, before the mod-$q$
book-keeping). Multiplying by $b$ shifts every remaining digit up one
place — $T[i+1]b^{m-1}+\dots+T[i+m-1]b^{1}$ — which is exactly $h_{i+1}$
minus its missing trailing term $T[i+m]b^0$. Adding $T[i+m]$ supplies that
term. $\blacksquare$

The constant $B = b^{m-1} \bmod q$ is the same for every window, so it is
computed **once**, via fast exponentiation ($O(\log m)$), and reused. Each
slide then costs $O(1)$: one subtraction, one multiplication, one addition,
all mod $q$. Building the first window's hash and the pattern's hash still
costs $O(m)$ each, so the whole scan is $O(n + m)$.

(Implementation note: `(a - c + q) % q` is used instead of `(a - c) % q`
throughout, because Java's `%` can return a negative result for a negative
left operand.)

## 4. The matching algorithm

At each position $i$: if $h_i = h(P)$, *verify* by comparing
$T[i \dots i+m-1]$ with $P$ character by character before reporting a
match. This is the **Las Vegas** form of the algorithm — it is always
exactly correct; a hash collision only wastes an $O(m)$ comparison, it
never produces a wrong answer. This is what `RabinKarp.search` implements.

A second, faster but riskier use of the same fingerprint — the **Monte
Carlo** form, from Rabin's original fingerprinting paper — skips
verification entirely and just compares $h(A) = h(B)$ to decide whether two
equal-length strings $A, B$ are equal. This is what `RabinKarp.probablyEqual`
implements; it is the right tool when you need to compare *many* long
strings cheaply (duplicate file detection, distributed replica checking,
etc.) and can tolerate a tiny, quantifiable error probability. Section 5
derives that probability.

## 5. Collision-probability analysis

This is the part that turns "seems to work" into a theorem.

**Setup.** Treat each length-$m$ string as an integer base $b$, so it lies
in $[0, b^m)$. Fix two *distinct* strings $A \ne B$ of length $m$ (for
`search`, $B$ ranges over the $n-m+1$ windows of $T$; for `probablyEqual`,
there is just one pair). Let

$$
d = |\,\text{int}(A) - \text{int}(B)\,|.
$$

Since $A \ne B$ and both lie in $[0, b^m)$, we have $0 < d < b^m$.

**Key fact.** $A$ and $B$ collide under modulus $q$ exactly when
$q \mid d$ (because $h(A) \equiv h(B) \pmod q \iff \mathrm{int}(A) \equiv
\mathrm{int}(B) \pmod q \iff q \mid d$). So a false collision happens
precisely when the randomly-chosen modulus happens to be one of $d$'s prime
factors.

**How many prime factors can $d$ have?** Since $d < b^m = 2^{m\log_2 b}$,
and the product of any $k$ *distinct* primes is at least $2^k$ (each factor
is $\ge 2$), a number less than $2^{m \log_2 b}$ can have at most
$m\log_2 b$ distinct prime factors. Call this bound $K = m\log_2 b$.

**Choosing $q$.** Pick $q$ uniformly at random from the primes in
$(2, M]$, for a bound $M$ we get to choose — crucially, chosen *before* and
*independently of* the strings being compared (this independence is what
makes the "Monte Carlo" guarantee hold even against an adversarially chosen
input; if $q$ were fixed in advance, an adversary could construct strings
that always collide under it). By the Prime Number Theorem there are
$\pi(M) \approx M/\ln M$ such primes. At most $K$ of them divide $d$, so

$$
\Pr[q \mid d] \;\le\; \frac{K}{\pi(M)} \;\approx\; \frac{m\log_2 b \cdot \ln M}{M}.
$$

**For pattern search** (`search`'s `probablyEqual`-free path only *uses*
the hash as a filter, so this bound is really about how many wasted
verifications occur, not correctness — correctness is already guaranteed
by the character check). Taking a union bound over the $n - m + 1$ windows
bounds the *expected number of spurious hash hits*:

$$
\mathbb{E}[\text{spurious hits}] \;\le\; (n-m+1)\cdot \frac{m\log_2 b \cdot \ln M}{M}
\;=\; O\!\left(\frac{nm\log b\,\log M}{M}\right).
$$

Choosing $M = \Theta(nm\log(nm))$ pushes this expectation below any
constant, and $q$ then needs only $O(\log(nm))$ bits — a handful more bits
than the input size itself, which is why in practice a single machine word
(or a modulus around $10^9$, as used in the code and tests below) is
already extremely safe for texts of realistic size.

**For `probablyEqual`** (single-pair Monte Carlo comparison, no
verification, no union bound needed): the false-positive probability for
one fixed pair $A \ne B$ is bounded directly by

$$
\Pr[\text{false ``equal''}] \;\le\; \frac{m\log_2 b}{\pi(M)}.
$$

E.g. with $b = 256$, $m = 1000$, and $M \approx 2^{31}$ (comfortably inside
a `long`), $\pi(M) \approx M/\ln M \approx 1.0\times 10^8$, giving a
collision bound on the order of $10^{-4}$ per comparison — and it shrinks
further, roughly linearly, as $M$ grows, since $\pi(M)$ grows almost
linearly with $M$ while $m \log_2 b$ stays fixed.

**Takeaway.** The randomness is in the *choice of $q$*, not in the input.
That is exactly why `withRandomPrimeModulus` draws $q$ uniformly from a
range of primes with `SecureRandom`, instead of hard-coding a fixed
modulus: a fixed modulus can be defeated by a specific adversarial input,
but a modulus unknown to the adversary and drawn independently of the data
cannot.

## 6. Complexity summary

| Phase | Time |
|---|---|
| Precompute pattern hash + $B = b^{m-1}\bmod q$ | $O(m + \log m)$ |
| Scan text, one $O(1)$ update per position | $O(n)$ |
| Verification of each hash hit | $O(m)$ per hit (expected $O(1)$ hits total, see §5) |
| **Total (Las Vegas `search`)** | **$O(n+m)$ expected, $O(nm)$ worst case** |
| **`probablyEqual`** | $O(m)$, with the error bound of §5 |

## 7. Implementation and tests

`RabinKarp.java` implements exactly the algorithm derived above:

- `fingerprint(s)` — Horner's-rule hash, §2.
- `search(text, pattern)` — the $O(1)$-per-step rolling scan of §3–4,
  with mandatory verification (Las Vegas).
- `probablyEqual(a, b)` — the unverified Monte-Carlo comparison of §4,
  whose error is bounded in §5.
- `withRandomPrimeModulus(base, bound)` — draws $q$ the way §5 requires:
  uniformly from the primes below `bound`, independent of any input.

The file is self-contained and runs its own test suite:

```
java RabinKarp.java
```

Tests included (all passing):

1. **Basic match** — the classic `"ABABDABACDABABCABAB"` / `"ABABCABAB"`
   example (expected hit at index 10).
2. **No match.**
3. **Overlapping matches** (`"AAAAA"` / `"AA"`).
4. **Pattern longer than text.**
5. **Pattern equals text.**
6. **Empty pattern** (defined here to report no matches).
7. **Forced collisions** — modulus deliberately set to `7` so real hash
   collisions occur, checked against a brute-force reference to confirm
   the Las Vegas verification step keeps the answer exact regardless.
8. **Randomized stress test** — 40 trials of random text/pattern pairs
   (some guaranteed to occur in the text, some not) with a randomly drawn
   prime modulus, each checked against brute force.
9. **Monte-Carlo fingerprint comparison** sanity checks for
   `probablyEqual`.

## 8. References

- R. M. Karp, M. O. Rabin, "Efficient randomized pattern-matching
  algorithms," *IBM Journal of Research and Development*, 31(2), 1987.
- cp-algorithms.com, "String Hashing" / "Rabin–Karp for String Matching."
