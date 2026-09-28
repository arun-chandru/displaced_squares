#!/usr/bin/env python3
"""Exact displaced-square-hole scalar tail: q>=4/5, x>=32.

Python 3.8+ and NumPy. Run:
    python -B -O verify_thin_tail.py [--output thin_tail_results.json]

The manuscript's midpoint and integer-sum estimates reduce the problem to
the closed rectangle 32<=x<=256, 0<=y<=8. Every acceptance decision is exact;
floating square roots only propose integers validated by squared inequalities.
All checks remain active under Python -O. No input files or network. The
default writes no files. The midpoint bound is independently certified by
the companion rectangular verifier; its analytic proof is in the manuscript.
"""

import argparse
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json
import numpy as np


D = 2048
H = 2**20
PI_LO, PI_HI = F(103993, 33102), F(104348, 33215)
CAP = 256*D


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def upward_roots(n):
    roots = np.sqrt(n).astype(np.int64)
    roots -= roots*roots > n
    roots += (roots+1)*(roots+1) <= n
    require(bool(np.all((roots*roots <= n) & ((roots+1)*(roots+1) > n))),
            'Every integer-root proposal is validated')
    return roots+(roots*roots < n)


def sum_enclosures():
    # F(z)=sum_(j>=0) sqrt(1-((j+1/2)/z)^2); G uses positive integer frequencies.
    midpoint = np.zeros(CAP+1, dtype=np.int64)
    integer = np.zeros(CAP+1, dtype=np.int64)
    require(4*CAP*CAP < 2**63 and 256*(2*CAP+1) < 2**63,
            'Radicand and unscaled sum range')
    for j in range(1, 513):
        first = (D*j+1)//2
        k = np.arange(first, CAP+1, dtype=np.int64)
        radicands = 4*k*k-D*D*j*j
        upper = upward_roots(radicands)
        if j % 2:
            midpoint[first:] += upper
        else:
            integer[first:] += upper
    require(int(max(midpoint.max(), integer.max()))*H+2*CAP < 2**63,
            'Fixed-point sum-product range')
    k = np.arange(1, CAP+1, dtype=np.int64)
    midpoint[1:] = (midpoint[1:]*H+2*k-1)//(2*k)
    integer[1:] = (integer[1:]*H+2*k-1)//(2*k)
    require(midpoint[0] == integer[0] == 0, 'Zero-width scalar convention')
    # Independent Python-isqrt checks at activation points and tail joins.
    for index in (0, D//2, D, 8*D, 32*D, CAP):
        exact_upper = [0, 0]
        for j in range(1, 513):
            rad = 4*index*index-D*D*j*j
            if rad >= 0:
                root = isqrt(rad)
                exact_upper[j % 2] += root+(root*root < rad)
        if index:
            expected_g = (exact_upper[0]*H+2*index-1)//(2*index)
            expected_f = (exact_upper[1]*H+2*index-1)//(2*index)
        else:
            expected_f = expected_g = 0
        require(int(midpoint[index]) == expected_f and int(integer[index]) == expected_g,
                'Independent Python-isqrt sum at grid index %s' % index)
    return midpoint, integer


def analytic_checks():
    def atan_bounds(d):
        value = sum(((-1)**k*F(1, (2*k+1)*d**(2*k+1))
                     for k in range(12)), F(0))
        return value, value+F(1, 25*d**25)
    al, ah = atan_bounds(5)
    bl, bh = atan_bounds(239)
    require(PI_LO < 16*al-4*bh < 16*ah-4*bl < PI_HI, 'Machin pi enclosures')
    tangent = F(1, 5)
    for _ in range(2):
        tangent = 2*tangent/(1-tangent*tangent)
    require((tangent-F(1, 239))/(1+tangent/F(239)) == 1,
            'Machin tangent identity')
    require(0 < 4*al-bh and 4*ah-bl < 1, 'Machin angle branch')
    require(F(2, 27) < F(3, 11)**2, 'C0 upper enclosure')
    for n, bound in ((8, F(71, 200)), (24, F(41, 200)), (32, F(177, 1000))):
        require(n*bound*bound > 1, 'Reciprocal-root upper enclosure')
    interior = F(3, 11)*(F(71, 200)+F(41, 200)+F(4, 5)*F(177, 1000))
    require(interior == F(2631, 13750) < F(1, 5), 'Interior analytic region')
    require(F(10, 7)**2 > 2, 'Upper enclosure for sqrt(2)')
    margin = F(1, 25)-F(3, 11)*(F(10, 7)+F(4, 5))/16
    require(margin == F(31, 15400) > 0, 'Displacement-uniform analytic tail')
    return dict(pi_lower=str(PI_LO), pi_upper=str(PI_HI),
                interior_x_onset=32, interior_y_minimum=8,
                interior_upper=str(interior), interior_target='1/5',
                interior_margin=str(F(1, 5)-interior),
                arbitrary_displacement_tail_onset=256,
                arbitrary_displacement_margin=str(margin))


def certify_finite_rectangle():
    f, g = sum_enclosures()
    boxes = np.array([(x*D, (x+1)*D, y*D, (y+1)*D)
                      for x in range(32, 256) for y in range(8)], dtype=np.int64)
    weight = 9*PI_LO*H/(4*D)
    require(weight.numerator*CAP < 2**63, 'Weyl-product range')
    require(H*256*D+CAP < 2**63, 'Floor-product range')
    require(14*int(max(f.max(), g.max()))+H*8+weight.numerator*CAP < 2**63,
            'Conservative signed box-comparison range')
    levels = []
    accepted = covered = tested = 0
    least_slack = None
    while len(boxes):
        l, u, v, w = boxes.T
        floor_term = (H*(u//D)*D+l-1)//l
        lower_weyl = (weight.numerator*l)//weight.denominator
        bound = 5*f[w]+5*f[u-v]+4*g[u]+floor_term-lower_weyl
        good = bound <= 0
        tested += len(boxes)
        accepted += int(good.sum())
        covered += int(((u[good]-l[good])*(w[good]-v[good])).sum())
        if np.any(good):
            slack = int((-bound[good]).min())
            least_slack = slack if least_slack is None else min(least_slack, slack)
        bad = boxes[~good]
        levels.append(dict(level=len(levels), tested=len(boxes),
                           accepted=int(good.sum()), remaining=len(bad)))
        if not len(bad):
            break
        require(len(levels) < 12, 'Subdivision cap: unresolved scalar boxes')
        l, u, v, w = bad.T
        require(bool(np.all((u-l) > 1) and np.all((w-v) > 1)),
                'Finest grid: unresolved scalar boxes')
        m, n = (l+u)//2, (v+w)//2
        boxes = np.concatenate((np.stack((l, m, v, n), axis=1),
                                np.stack((l, m, n, w), axis=1),
                                np.stack((m, u, v, n), axis=1),
                                np.stack((m, u, n, w), axis=1)))
    require(covered == 224*8*D*D, 'Exact closed-rectangle covering area')
    require(tested == 1792+4*(tested-accepted), 'Full quadtree accounting')
    require((tested, accepted, len(levels)-1, least_slack)
            == (1_153_748, 865_759, 8, 0), 'Deterministic scalar checkpoints')
    require(levels[-1]['remaining'] == 0, 'No unresolved scalar box')
    return dict(x_range=[32, 256], y_range=[0, 8],
                grid_denominator=D, fixed_point_scale=H,
                levels=levels, tested=tested, accepted_leaves=accepted,
                maximum_depth=len(levels)-1, unresolved_boxes=0,
                covered_integer_area=covered, least_integer_slack=least_slack,
                independent_integer_root_sum_checks=6)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Optional explicit JSON destination')
    args = parser.parse_args()
    result = dict(certificate_version=1, certificate='displaced_square_thin_tail',
                  arithmetic='Exact integers/rationals; validated floating-root proposals',
                  q_range='[4/5,1)', x_range='[32,infinity)',
                  analytic=analytic_checks(), finite=certify_finite_rectangle(),
                  all_checks_passed=True)
    encoded = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output is not None:
        args.output.write_text(encoded, encoding='utf-8')
    print(encoded, end='')


if __name__ == '__main__':
    main()
