"""Verify the two finite transport certificates for the PL+ theorem.

Requires Python 3 and NumPy.  Run: python -B -O verify_transport.py
The default invocation reads only wide_certificate.json and
thin_certificate.json beside this file, then prints a JSON summary.
An optional --output PATH writes that summary to PATH.

This verifies the finite-index part, conditional on the manuscript's actual
operator gap/transport lemmas and two analytic tails.  It does not use a
generator, stored PASS claims, stored replay results, or approximate spectra.
"""

import argparse
from fractions import Fraction
import json
from math import isqrt
from pathlib import Path

import numpy as np


CAP = 250000
SCALE = 2**40
PI_LO = Fraction(103993, 33102)
PI_HI = Fraction(104348, 33215)
INT64_LIMIT = 2**63
HERE = Path(__file__).resolve().parent


def require(condition, message):
    """Proof obligations remain active with Python -O."""
    if not condition:
        raise RuntimeError(message)


def integer(value, name):
    require(type(value) is int, name + ': expected an integer')
    return value


def floor_fraction(value):
    return value.numerator // value.denominator


def certify_pi():
    """Exact alternating series, using pi = 16 atan(1/5)-4 atan(1/239)."""
    def enclosure(x):
        terms = 12
        partial = sum(((-1)**j * x**(2*j+1) / (2*j+1)
                       for j in range(terms)), Fraction(0))
        following = (-1)**terms * x**(2*terms+1) / (2*terms+1)
        return min(partial, partial+following), max(partial, partial+following)

    low5, high5 = enclosure(Fraction(1, 5))
    low239, high239 = enclosure(Fraction(1, 239))
    require(PI_LO < 16*low5 - 4*high239, 'Certified strict lower pi bound')
    require(16*high5 - 4*low239 < PI_HI, 'Certified strict upper pi bound')


def disk_counts():
    """Count every ordered positive pair, with no symmetry reduction."""
    counts = np.zeros(CAP+1, dtype=np.int64)
    for m in range(1, isqrt(CAP)+1):
        ns = np.arange(1, isqrt(CAP-m*m)+1, dtype=np.int64)
        # For fixed m, all destinations are distinct; += loses no repetitions.
        counts[m*m + ns*ns] += 1
    np.cumsum(counts, out=counts)
    require(0 <= int(counts[-1]) <= CAP, 'Disk-count range')
    for n in (0, 1, 2, 40320, 40401, 65536, CAP):
        direct = sum(isqrt(n-m*m) for m in range(1, isqrt(n)+1))
        require(int(counts[n]) == direct, 'Independent radial count at '+str(n))
    return counts


def transport(q, r):
    if q >= r:
        return Fraction(1)
    return q*(1-r)**3 / (r*(1-q)**3)


def coefficient(left, right, reference):
    endpoint_minimum = min((1-q*q)*transport(q, reference)
                           for q in (left, right))
    return floor_fraction(PI_LO * SCALE * endpoint_minimum / 4)


def verify_certificate(path, band, counts):
    require(band in ('wide', 'thin'), 'Known band')
    result = json.loads(Path(path).read_text(encoding='utf-8-sig'))
    if band == 'wide':
        denominator = 64000
        first, last = Fraction(12, 25), Fraction(4, 5)
        expected_intervals, expected_repairs = 582, 130
        require(result.get('onset_radius') == 201, 'Wide certificate onset')
    else:
        denominator = 192000
        first, last = Fraction(4, 5), Fraction(11, 12)
        expected_intervals, expected_repairs = 739, 81

    require(integer(result['count_cap'], 'Count cap') == CAP, 'Fixed count cap')
    require(integer(result['scale'], 'Scale') == SCALE, 'Fixed integer scale')
    require(integer(result['reference_denominator'], 'Reference denominator')
            == denominator, 'Fixed reference denominator')
    require(Fraction(result['pi_lower']) == PI_LO
            and Fraction(result['pi_upper']) == PI_HI, 'Canonical pi bounds')
    require(SCALE*CAP < INT64_LIMIT
            and denominator**2*CAP < INT64_LIMIT, 'All vector products fit int64')

    records = result['records']
    require(isinstance(records, list) and len(records) > 0, 'Nonempty records')
    prepared = []
    for record in records:
        interval = record['q_interval']
        require(isinstance(interval, list) and len(interval) == 2, 'Two interval endpoints')
        left, right = map(Fraction, interval)
        require(first <= left < right <= last, 'Valid closed interval')
        prepared.append((left, right, record))
    prepared.sort(key=lambda row: row[0])
    require(prepared[0][0] == first and prepared[-1][1] == last, 'Exact band endpoints')
    for previous, following in zip(prepared, prepared[1:]):
        require(previous[1] == following[0], 'Consecutive closed intervals meet exactly')

    tau = np.arange(1, CAP+1, dtype=np.int64)

    def frontier(numerator, kmax):
        inner = numerator*numerator*tau // (denominator*denominator)
        indices = counts[tau-1] - counts[inner] + 1
        require(bool(np.all(indices >= 1)), 'Positive spectral-gap indices')
        available = indices <= kmax
        lower = np.zeros(kmax+1, dtype=np.int64)
        np.maximum.at(lower, indices[available], tau[available])
        np.maximum.accumulate(lower, out=lower)
        return lower

    total_pairs = 0
    total_repairs = 0
    maximum_kmax = 0
    for left, right, record in prepared:
        energy_cutoff = (Fraction(201**2)/(left*left) if band == 'wide'
                         else Fraction(32**2)/(1-right)**2)
        kmax = floor_fraction(PI_HI*(1-left*left)*energy_cutoff/4)
        require(integer(record['kmax'], 'Index cutoff') == kmax
                and 1 <= kmax < CAP, 'Exact conservative finite-index cutoff')
        references = record['default_references']
        coefficients = record['default_scaled_coefficients']
        require(isinstance(references, list) and isinstance(coefficients, list)
                and len(references) == len(coefficients) == 2, 'Two paired default references')
        ks = np.arange(1, kmax+1, dtype=np.int64)
        passed = np.zeros(kmax, dtype=bool)
        for numerator, saved_coefficient in zip(references, coefficients):
            numerator = integer(numerator, 'Reference numerator')
            require(0 < numerator < denominator, 'Strictly interior reference size')
            reference = Fraction(numerator, denominator)
            mu = coefficient(left, right, reference)
            require(integer(saved_coefficient, 'Default coefficient') == mu
                    and 0 <= mu < SCALE, 'Exact conservative default coefficient')
            passed |= mu*frontier(numerator, kmax)[1:] >= SCALE*ks

        failures = ks[~passed].tolist()
        saved_failures = record['default_failed_indices']
        require(isinstance(saved_failures, list), 'Default failures are a list')
        for k in saved_failures:
            integer(k, 'Failed index')
        require(saved_failures == failures, 'Every default outcome reconstructed')
        repairs = record['repairs']
        require(isinstance(repairs, list), 'Repairs are a list')
        repaired = []
        for repair in repairs:
            k = integer(repair['k'], 'Repair index')
            numerator = integer(repair['reference'], 'Repair reference')
            n = integer(repair['tau'], 'Repair energy')
            require(k in failures and 0 < numerator < denominator
                    and 1 <= n <= CAP, 'Repair scope')
            inner = numerator*numerator*n // (denominator*denominator)
            gap_index = int(counts[n-1]) - int(counts[inner]) + 1
            require(1 <= gap_index <= k
                    and integer(repair['gap_index'], 'Saved gap index') == gap_index,
                    'Strict-outer inclusive-inner spectral gap')
            mu = coefficient(left, right, Fraction(numerator, denominator))
            require(integer(repair['scaled_coefficient'], 'Repair coefficient') == mu
                    and 0 <= mu < SCALE and mu*n >= SCALE*k, 'Repair reaches target')
            repaired.append(k)
        require(sorted(repaired) == failures, 'Every default failure repaired exactly once')
        total_pairs += kmax
        total_repairs += len(repaired)
        maximum_kmax = max(maximum_kmax, kmax)

    # These fixed totals identify the released data; coverage was checked above.
    require(len(prepared) == expected_intervals, 'Released interval total')
    require(total_repairs == expected_repairs, 'Released repair total')
    return {'interval': [str(first), str(last)], 'reference_denominator': denominator,
            'intervals': len(prepared), 'interval_index_pairs': total_pairs,
            'repairs': total_repairs, 'maximum_kmax': maximum_kmax}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Optionally save the JSON summary')
    args = parser.parse_args()
    certify_pi()
    counts = disk_counts()
    summary = {
        'status': 'PASS_FINITE_TRANSPORT_CERTIFICATES',
        'optimized_python': not __debug__,
        'count_cap': CAP,
        'independent_last_count': int(counts[-1]),
        'pi_lower': str(PI_LO), 'pi_upper': str(PI_HI),
        'wide': verify_certificate(HERE/'wide_certificate.json', 'wide', counts),
        'thin': verify_certificate(HERE/'thin_certificate.json', 'thin', counts),
        'mathematical_dependencies': [
            'Actual Dirichlet gap and transport lemmas in the manuscript.',
            'Wide: 8u/9 <= d(u) <= 10u/9 for all real u >= 201; this supplies the tail for t >= 201/q.',
            'Thin: mixed scalar deficit D_q(x,a) <= 0 for q >= 4/5, x >= 32, a in [0,1]; this supplies the tail for t >= 32/(1-q).'
        ]
    }
    text = json.dumps(summary, indent=2)+'\n'
    if args.output is not None:
        args.output.write_text(text, encoding='utf-8')
    print(text, end='')


if __name__ == '__main__':
    main()
