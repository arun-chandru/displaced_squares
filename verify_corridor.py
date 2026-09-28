#!/usr/bin/env python3
"""Exact disk census proving 8u/9 <= d(u) <= 10u/9 for every u>=201.

Python 3.8+ and NumPy. Run:
    python -B -O verify_corridor.py [--output corridor_results.json]

The manuscript's explicit disk-discrepancy estimate supplies the analytic
tail from radius 15625. This program checks its constants and every finite
squared-radius band. Float square roots only propose integers, which are
validated by exact inequalities before use. Acceptance is exact and remains
active under Python -O. No input files or network. Default: no file writes.
"""

import argparse
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json
import numpy as np


RADIUS = 15625
CAP = RADIUS**2
SEGMENT = 4_000_000
ONSET = 201
PI_LO, PI_HI = F(103993, 33102), F(104348, 33215)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def atan_bounds(d):
    value = sum(((-1)**k*F(1, (2*k+1)*d**(2*k+1))
                 for k in range(12)), F(0))
    return value, value+F(1, 25*d**25)


def certify_constants():
    al, ah = atan_bounds(5)
    bl, bh = atan_bounds(239)
    require(PI_LO < 16*al-4*bh < 16*ah-4*bl < PI_HI, 'Machin pi enclosures')
    tangent = F(1, 5)
    for _ in range(2):
        tangent = 2*tangent/(1-tangent*tangent)
    require((tangent-F(1, 239))/(1+tangent/F(239)) == 1,
            'Machin tangent identity')
    require(0 < 4*al-bh and 4*ah-bl < 1, 'Machin angle branch')
    require(RADIUS == 25**3 and CAP == 244_140_625, 'Analytic tail onset')
    lower = 1-F(3, 4*RADIUS)-F(11, 100)
    upper = 1+F(1, 4*RADIUS)+F(11, 100)
    require(lower == F(27811, 31250) > F(8, 9), 'Analytic lower tail')
    require(upper == F(17344, 15625) < F(10, 9), 'Analytic upper tail')
    # C_N <= floor(sqrt(N))^2 <= N, so every cumulative uint32 count fits.
    require(CAP < 2**32, 'Unsigned cumulative count range')
    require(9*PI_HI.numerator*(CAP+1)+36*PI_HI.denominator*CAP < 2**63,
            'Signed comparison-product range')
    require(40*PI_HI.denominator*(RADIUS+1) < 2**63,
            'Signed endpoint-product range')
    return dict(pi_lower=str(PI_LO), pi_upper=str(PI_HI),
                analytic_radius_onset=RADIUS, analytic_lower=str(lower),
                analytic_upper=str(upper))


def integer_roots(n):
    roots = np.sqrt(n).astype(np.int64)
    roots -= roots*roots > n
    roots += (roots+1)*(roots+1) <= n
    require(bool(np.all((roots*roots <= n) & ((roots+1)*(roots+1) > n))),
            'Every integer-root proposal is validated')
    return roots


def row_count(n):
    return sum(isqrt(n-a*a) for a in range(1, isqrt(n)+1))


def certify_corridor():
    carry = processed = 0
    bad_lower, bad_upper = [], []
    checkpoints = {ONSET**2: None, 1_000_000: None, CAP: None}
    for begin in range(0, CAP+1, SEGMENT):
        end = min(CAP+1, begin+SEGMENT)
        counts = np.zeros(end-begin, dtype=np.uint32)
        for m in range(1, isqrt(end-1)+1):
            mm = m*m
            rem = max(0, begin-mm)
            first = isqrt(rem)
            first = max(1, first+(first*first < rem))
            last = min(m-1, isqrt(end-1-mm))
            if first <= last:
                n = np.arange(first, last+1, dtype=np.int64)
                # n^2 is injective in this row; advanced indexing has no repeats.
                counts[mm+n*n-begin] += 2
            if begin <= 2*mm < end:
                counts[2*mm-begin] += 1
        np.cumsum(counts, dtype=np.uint32, out=counts)
        counts += np.uint32(carry)
        carry = int(counts[-1])
        require(carry <= end-1, 'Cumulative lattice-count bound')
        for radius_squared in checkpoints:
            if begin <= radius_squared < end:
                checkpoints[radius_squared] = int(counts[radius_squared-begin])
        for offset in range(0, end-begin, 65536):
            finish = min(offset+65536, end-begin)
            n = np.arange(begin+offset, begin+finish, dtype=np.int64)
            count = counts[offset:finish].astype(np.int64)
            roots = integer_roots(n)
            upper_roots = roots+(roots*roots < n)
            next_roots = integer_roots(n+1)
            lower = PI_LO.numerator*n-4*PI_LO.denominator*count
            upper = PI_HI.numerator*(n+1)-4*PI_HI.denominator*count
            failed_lower = 9*lower < 32*PI_LO.denominator*upper_roots
            failed_upper = 9*upper > 40*PI_HI.denominator*next_roots
            bad_lower.extend(n[failed_lower].tolist())
            bad_upper.extend(n[failed_upper].tolist())
        processed += end-begin
    require(processed == CAP+1, 'Every squared-radius band processed')
    last_bad = max(bad_lower+bad_upper, default=0)
    derived_onset = isqrt(last_bad)+1
    require(derived_onset <= ONSET, 'Claimed corridor onset 201')
    require((len(bad_lower), max(bad_lower), len(bad_upper), max(bad_upper))
            == (498, 17226, 481, 40320), 'Complete failure-list checkpoints')
    require(carry == 191_731_890, 'Terminal count checkpoint')
    for n, observed in checkpoints.items():
        require(observed == row_count(n), 'Independent row count at n=%s' % n)
    return dict(radius_onset=ONSET, derived_onset=derived_onset,
                cap_squared_radius=CAP, processed_squared_bands=processed,
                segment_bytes=SEGMENT*4, count_at_cap=carry,
                lower_bad_count=len(bad_lower), lower_bad_max=max(bad_lower),
                upper_bad_count=len(bad_upper), upper_bad_max=max(bad_upper),
                lower_bad=bad_lower, upper_bad=bad_upper,
                independent_row_counts={str(n): c for n, c in checkpoints.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Optional explicit JSON destination')
    args = parser.parse_args()
    result = dict(certificate_version=1, certificate='disk_deficit_corridor',
                  arithmetic='Exact integers/rationals; validated floating-root proposals',
                  constants=certify_constants(), corridor=certify_corridor(),
                  all_checks_passed=True)
    encoded = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output is not None:
        args.output.write_text(encoded, encoding='utf-8')
    print(encoded, end='')


if __name__ == '__main__':
    main()
