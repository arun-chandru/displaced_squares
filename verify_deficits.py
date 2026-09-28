#!/usr/bin/env python3
"""Exact square-deficit and multiple-square-hole certificates for PL+.

Python 3.8+; standard library only. Run:
    python -B -O verify_deficits.py [--output deficit_results.json]

The manuscript proves the analytic disk-discrepancy estimate and interval
reductions. This program checks their explicit constants and every finite
scalar obligation. It does not approximate perforated-domain eigenvalues.
All checks are unconditional, including under Python -O. No input files,
network, or floating-point arithmetic are used. The default writes no files.
"""

import argparse
from array import array
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json


CAP = 1_000_000
PI_LO, PI_HI = F(333, 106), F(355, 113)
SIDE_LO, SIDE_HI = F(97, 200), F(243, 500)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def atan_enclosure(q, terms):
    require(q > 1 and terms > 0, 'Invalid arctangent parameters')
    value = sum((F((-1)**k, (2*k+1)*q**(2*k+1))
                 for k in range(terms)), F(0))
    following = value + F((-1)**terms, (2*terms+1)*q**(2*terms+1))
    return min(value, following), max(value, following)


def certify_constants():
    al, ah = atan_enclosure(5, 12)
    bl, bh = atan_enclosure(239, 4)
    machin_lo, machin_hi = 16*al-4*bh, 16*ah-4*bl
    require(PI_LO < machin_lo < machin_hi < PI_HI, 'Rational pi enclosure')
    tangent = F(1, 5)
    for _ in range(2):
        tangent = 2*tangent/(1-tangent*tangent)
    require((tangent-F(1, 239))/(1+tangent/F(239)) == 1,
            'Machin tangent identity')
    require(0 < 4*al-bh and 4*ah-bl < 1, 'Machin angle branch')

    s2_lo, s2_hi = 1-F(12, 5)/PI_LO, 1-F(12, 5)/PI_HI
    require(SIDE_LO**2 < s2_lo < s2_hi < SIDE_HI**2,
            'Sharp side threshold enclosure')
    transitions = []
    for numerator, left, right in ((2, 8, 9), (5, 21, 22), (18, 76, 77)):
        require(left < F(numerator)/s2_hi <= F(numerator)/s2_lo < right,
                'Envelope transition enclosure')
        transitions.append(dict(numerator=numerator, strict_lower=left,
                                strict_upper=right))
    require(PI_HI**2*F(2, 16) < F(1309, 1000)**2 < F(4, 3)**2,
            'Increasing envelope slopes')
    require(0 < s2_hi < F(2, 5), 'Sharpness witness interval')

    require(F(16, 45) < F(9, 25) and F(18, 25) < 1,
            'Bessel bound on [1,2]')
    require(F(32, 13)/PI_LO < 1, 'Bessel energy constant')
    require(F(157, 50) < PI_LO and PI_HI < F(22, 7), 'Coarser pi bounds')
    require(F(1, 16) < F(2, 5)**3, 'Two to the power -4/3')
    require(F(19, 12)**3 < 4, 'Two to the power 2/3')
    require(F(12, 5) < F(31, 20)**2, 'Square-root majorant')
    discrepancy = F(121, 35)+F(22320, 2983)
    require(discrepancy == F(1142143, 104405) < 11,
            'Analytic discrepancy coefficient')
    lower1000 = 1-F(3, 4000)-F(11, 40)
    upper1000 = 1+F(1, 4000)+F(11, 40)
    require(CAP == 1000**2 and 1000 == 10**3,
            'Finite range joins the analytic radius-1000 tail')
    require(lower1000 == F(2897, 4000) > F(81, 125),
            'Analytic lower slope at radius 1000')
    require(upper1000 == F(5101, 4000) < F(4, 3),
            'Analytic upper slope at radius 1000')
    return dict(pi_lower=str(PI_LO), pi_upper=str(PI_HI),
                side_threshold='sqrt(1-12/(5*pi))',
                side_lower=str(SIDE_LO), side_upper=str(SIDE_HI),
                squared_transition_intervals=transitions,
                discrepancy_constant_upper=str(discrepancy),
                discrepancy_radius_onset=2,
                deficit_tail_radius_onset=1000,
                radius1000_lower_slope=str(lower1000),
                radius1000_upper_slope=str(upper1000))


def build_square_counts():
    counts = array('I', [0])*(CAP+1)
    require(counts.itemsize >= 4 and CAP < 2**32,
            'Unsigned count storage range')
    # Each positive ordered pair is counted once, including multiplicity.
    for a in range(1, isqrt(CAP//2)+1):
        aa = a*a
        for b in range(a, isqrt(CAP-aa)+1):
            counts[aa+b*b] += 1 if a == b else 2
    for n in range(1, CAP+1):
        counts[n] += counts[n-1]
    require(counts[CAP] == 784_387, 'Count at squared radius 1,000,000')
    return counts


def row_count(n):
    return sum(isqrt(n-a*a) for a in range(1, isqrt(n)+1))


def certify_deficits(counts):
    exceptions = []
    for n in range(1, CAP+1):
        # Left limits for upper bounds; inclusive values for lower bounds.
        upper = 355*n-452*counts[n-1]
        require(upper <= 0 or 9*upper*upper <= 1808**2*n,
                'Global upper deficit at n=%s' % n)
        if n <= 18:
            require(upper <= 0 or (1000*upper)**2 <= (1309*452)**2*n,
                    'Prefix upper deficit at n=%s' % n)
        lower = 333*n-424*counts[n]
        if 21 <= n <= 76:
            if lower < 0 or (500000*lower)**2 < (318087*424)**2*n:
                exceptions.append(n)
        if 76 <= n:
            require(lower >= 0 and (125*lower)**2 >= (81*424)**2*n,
                    'High lower deficit at n=%s' % n)
    require(exceptions == [53], 'Complete exceptional-band list')
    for n in range(9):
        require(3*n >= 5*counts[n], 'First quadratic envelope branch')
    for n in range(8, 22):
        lower = 333*n-424*counts[n]
        require(3*n >= 5*(counts[n]-1), 'Second quadratic envelope branch')
        require(lower >= 0 and (500*lower)**2 >= 2*n*(243*333)**2,
                'Intermediate linear envelope branch')
    independent_radii = list(range(77))+[CAP-1, CAP]
    for n in independent_radii:
        require(counts[n] == row_count(n), 'Independent row count at n=%s' % n)
    return dict(cap_squared_radius=CAP, count_at_cap=counts[CAP],
                global_upper_left_limit_range=[1, CAP],
                prefix_upper_left_limit_range=[1, 18],
                middle_lower_range=[21, 76], exceptional_bands=exceptions,
                high_lower_range=[76, CAP],
                global_upper_left_limit_failures=[],
                prefix_upper_left_limit_failures=[], high_lower_failures=[],
                low_branch_failures=[],
                independent_row_count_checks=len(independent_radii),
                low_counts_n_0_through_76=list(counts[:77]))


def certify_partition(counts):
    budget, ceiling = F(18, 5), F(9, 2)
    require(SIDE_HI**2*54 < budget**2 < 13, 'Exceptional-band side budget')
    require(F(3, 2)**2 > 2, 'Upper bound for sqrt(2)')
    small = PI_HI*F(3, 2)*budget/4
    require(small == F(1917, 452) < ceiling, 'All-small partition')
    endpoints = []
    for n, root_lo, count, expected in (
        (2, F(7, 5), 1, F(2488, 565)),
        (5, F(223, 100), 1, F(49973, 11300)),
        (8, F(14, 5), 3, F(442, 113)),
        (10, F(79, 25), 4, F(11349, 2825)),
    ):
        require(0 < root_lo and root_lo**2 < n, 'Positive radical lower bound')
        actual = counts[n] if n == 2 else counts[n-1]
        require(actual == count, 'Inclusive/left-limit partition count')
        bound = PI_HI*(2*n+budget**2-2*budget*root_lo)/4-count
        require(bound == expected < ceiling, 'Partition endpoint bound')
        endpoints.append(dict(squared_radius=n, radical_lower=str(root_lo),
                              count=count, upper=str(bound)))
    require(counts[12] == 6 and counts[53] == 37, 'Closing partition counts')
    last, outer = PI_HI*budget**2/4-6, 53*PI_LO/4-counts[53]
    require(last == F(2361, 565) < ceiling < outer == F(37, 8),
            'Exceptional-band closing inequalities')
    return dict(budget=str(budget), deficit_ceiling=str(ceiling),
                all_small_upper=str(small), endpoint_bounds=endpoints,
                final_endpoint_upper=str(last), outer_53_lower=str(outer))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Optional explicit JSON destination')
    args = parser.parse_args()
    constants = certify_constants()
    counts = build_square_counts()
    result = dict(certificate_version=1, certificate='square_deficits',
                  arithmetic='Python integers, integer isqrt, fractions.Fraction',
                  constants=constants, deficits=certify_deficits(counts),
                  partition=certify_partition(counts), all_checks_passed=True)
    encoded = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output is not None:
        args.output.write_text(encoded, encoding='utf-8')
    print(encoded, end='')


if __name__ == '__main__':
    main()
