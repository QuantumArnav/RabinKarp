import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;

/**
 * Rabin-Karp / Karp-Rabin fingerprinting.
 *
 * A window of m characters starting at text position i is treated as an
 * m-digit number in base `base`:
 *
 *     h_i = ( T[i]*base^(m-1) + T[i+1]*base^(m-2) + ... + T[i+m-1]*base^0 )  mod q
 *
 * Sliding the window one place to the right can be done in O(1):
 *
 *     h_{i+1} = ( (h_i - T[i]*base^(m-1)) * base + T[i+m] )  mod q
 *
 * Two independent uses of the same fingerprint are implemented here:
 *
 *   1. search()        -- "Las Vegas" substring search: a hash match is
 *                          always verified character-by-character before
 *                          being reported, so the answer is exactly correct.
 *                          A collision only costs a wasted O(m) comparison.
 *
 *   2. probablyEqual()  -- "Monte Carlo" fingerprint comparison in the
 *                          sense of Rabin (1981): two strings are declared
 *                          equal iff their fingerprints agree, with NO
 *                          verification. This can be wrong, but only with
 *                          the probability bounded in the accompanying
 *                          report (roughly n*m*log(base) / q for an
 *                          adversarial input and a modulus q chosen
 *                          uniformly at random and independently of the
 *                          input).
 *
 * Run this file directly to execute the self-test suite:
 *     java RabinKarp.java
 */
public class RabinKarp {

    private final long base;
    private final long mod;

    public RabinKarp(long base, long mod) {
        if (mod <= 1) throw new IllegalArgumentException("modulus must be > 1");
        if (base <= 0) throw new IllegalArgumentException("base must be > 0");
        this.base = base;
        this.mod = mod;
    }

    /** Builds an instance whose modulus is a prime chosen uniformly at
     *  random from (2, modulusBound), independently of any input string --
     *  this is the randomization the collision-probability analysis
     *  depends on. */
    public static RabinKarp withRandomPrimeModulus(long base, long modulusBound) {
        return new RabinKarp(base, randomPrimeBelow(modulusBound, new SecureRandom()));
    }

    public long getBase() { return base; }
    public long getModulus() { return mod; }

    /** h(s) mod q, computed in O(|s|) via Horner's rule. */
    public long fingerprint(String s) {
        long h = 0;
        for (int k = 0; k < s.length(); k++) {
            h = (h * base + s.charAt(k)) % mod;
        }
        return h;
    }

    /** All start indices i such that text.substring(i, i+m) equals pattern.
     *  Always exact (verified), regardless of how collision-prone `mod` is. */
    public List<Integer> search(String text, String pattern) {
        List<Integer> hits = new ArrayList<>();
        int n = text.length(), m = pattern.length();
        if (m == 0 || m > n) return hits;

        long hPattern = fingerprint(pattern);
        long hWindow = fingerprint(text.substring(0, m));
        long highOrder = powmod(base, m - 1, mod); // base^(m-1) mod q, computed once

        int i = 0;
        while (true) {
            if (hWindow == hPattern && text.regionMatches(i, pattern, 0, m)) {
                hits.add(i);
            }
            int next = i + m;
            if (next >= n) break;
            long leading = (text.charAt(i) % mod) * highOrder % mod;
            hWindow = (hWindow - leading + mod) % mod;      // drop the outgoing digit
            hWindow = (hWindow * base + text.charAt(next)) % mod; // shift in the new one
            i++;
        }
        return hits;
    }

    /** Monte-Carlo equality check: same length and same fingerprint.
     *  No character-by-character verification is performed. */
    public boolean probablyEqual(String a, String b) {
        return a.length() == b.length() && fingerprint(a) == fingerprint(b);
    }

    // ---------------- internals ----------------

    private static long powmod(long b, long e, long m) {
        long result = 1 % m, base = b % m;
        while (e > 0) {
            if ((e & 1) == 1) result = (result * base) % m;
            base = (base * base) % m;
            e >>= 1;
        }
        return result;
    }

    private static long randomPrimeBelow(long bound, SecureRandom rnd) {
        if (bound < 5) throw new IllegalArgumentException("bound too small");
        while (true) {
            long c = 3 + (long) (rnd.nextDouble() * (bound - 3));
            if (c % 2 == 0) c++;
            if (BigInteger.valueOf(c).isProbablePrime(40)) return c;
        }
    }

    private static List<Integer> bruteForceSearch(String text, String pattern) {
        List<Integer> hits = new ArrayList<>();
        int n = text.length(), m = pattern.length();
        if (m == 0 || m > n) return hits;
        for (int i = 0; i <= n - m; i++) {
            if (text.regionMatches(i, pattern, 0, m)) hits.add(i);
        }
        return hits;
    }

    // ---------------- self-tests ----------------

    public static void main(String[] args) {
        List<Runnable> tests = List.of(
                RabinKarp::testBasicMatch,
                RabinKarp::testNoMatch,
                RabinKarp::testOverlappingMatches,
                RabinKarp::testPatternLongerThanText,
                RabinKarp::testPatternEqualsText,
                RabinKarp::testEmptyPattern,
                RabinKarp::testForcedCollisionsStillCorrect,
                RabinKarp::testRandomizedAgainstBruteForce,
                RabinKarp::testFingerprintEqualityMonteCarlo
        );

        int passed = 0, failed = 0;
        for (Runnable t : tests) {
            try {
                t.run();
                passed++;
            } catch (AssertionError e) {
                failed++;
                System.out.println("FAIL: " + e.getMessage());
            }
        }
        System.out.println(passed + " passed, " + failed + " failed");
        if (failed > 0) System.exit(1);
    }

    private static void check(boolean cond, String msg) {
        if (!cond) throw new AssertionError(msg);
    }

    private static void testBasicMatch() {
        RabinKarp rk = new RabinKarp(256, 1_000_000_007L);
        List<Integer> hits = rk.search("ABABDABACDABABCABAB", "ABABCABAB");
        check(hits.equals(List.of(10)), "expected match at index 10, got " + hits);
    }

    private static void testNoMatch() {
        RabinKarp rk = new RabinKarp(256, 1_000_000_007L);
        check(rk.search("abcdefgh", "xyz").isEmpty(), "expected no matches");
    }

    private static void testOverlappingMatches() {
        RabinKarp rk = new RabinKarp(256, 1_000_000_007L);
        List<Integer> hits = rk.search("AAAAA", "AA");
        check(hits.equals(List.of(0, 1, 2, 3)), "expected overlapping matches, got " + hits);
    }

    private static void testPatternLongerThanText() {
        RabinKarp rk = new RabinKarp(256, 1_000_000_007L);
        check(rk.search("ab", "abc").isEmpty(), "pattern longer than text must give no matches");
    }

    private static void testPatternEqualsText() {
        RabinKarp rk = new RabinKarp(256, 1_000_000_007L);
        check(rk.search("hello", "hello").equals(List.of(0)), "pattern == text should match once at 0");
    }

    private static void testEmptyPattern() {
        RabinKarp rk = new RabinKarp(256, 1_000_000_007L);
        check(rk.search("hello", "").isEmpty(), "empty pattern is defined to return no matches here");
    }

    /** Deliberately use a tiny modulus so real hash collisions occur, and
     *  confirm the Las-Vegas verification step still gives an exactly
     *  correct answer -- this is the practical payoff of the
     *  collision-probability analysis in the report. */
    private static void testForcedCollisionsStillCorrect() {
        RabinKarp rk = new RabinKarp(256, 7); // tiny modulus -> frequent collisions
        String text = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBB";
        String pattern = "AAAB";
        check(rk.search(text, pattern).equals(bruteForceSearch(text, pattern)),
                "search must match brute force even with a collision-prone modulus");
    }

    private static void testRandomizedAgainstBruteForce() {
        SecureRandom rnd = new SecureRandom();
        for (int trial = 0; trial < 40; trial++) {
            int n = 100 + rnd.nextInt(300);
            int m = 1 + rnd.nextInt(Math.min(15, n));

            StringBuilder sb = new StringBuilder();
            for (int k = 0; k < n; k++) sb.append((char) ('a' + rnd.nextInt(4))); // small alphabet -> frequent near-matches
            String text = sb.toString();

            String pattern;
            if (trial % 2 == 0) {
                int start = rnd.nextInt(n - m + 1);
                pattern = text.substring(start, start + m); // guaranteed to occur at least once
            } else {
                StringBuilder pb = new StringBuilder();
                for (int k = 0; k < m; k++) pb.append((char) ('a' + rnd.nextInt(4)));
                pattern = pb.toString(); // may or may not occur
            }

            RabinKarp rk = withRandomPrimeModulus(256, 1_000_000);
            List<Integer> got = rk.search(text, pattern);
            List<Integer> want = bruteForceSearch(text, pattern);
            check(got.equals(want), "trial " + trial + ": mismatch for pattern \"" + pattern
                    + "\" got=" + got + " want=" + want);
        }
    }

    private static void testFingerprintEqualityMonteCarlo() {
        RabinKarp rk = new RabinKarp(256, 1_000_000_007L);
        check(rk.probablyEqual("same string", "same string"), "identical strings must compare equal");
        check(!rk.probablyEqual("same string", "different string!"), "different-length strings must differ");
        check(!rk.probablyEqual("abcdef", "abcdeg"), "distinct same-length strings should (almost always) differ");
    }
}