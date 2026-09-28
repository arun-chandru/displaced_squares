# Pólya inequalities for displaced square holes and perforated squares

**Author:** R. Arun Chandru  
**Date:** 28 September 2026  
Replacement for *Pólya inequalities for square frames and perforated squares* (27 September 2026).

## Abstract

We prove Pólya's inequality for every Dirichlet eigenvalue of a square with an arbitrary strictly interior parallel square removed. The hole may have any size and any position; neither centering nor a positive lower bound on the boundary gaps is required. The proof combines complementary square counts, endpoint-sensitive spectral gaps, and a displacement-uniform transport between the actual Dirichlet form domains. Mixed-boundary comparisons and explicit scalar estimates cover the thin-gap regime. We also prove the inequality for finitely many independently positioned and oriented square holes whose total side length is at most `sqrt(1 - 12/(5*pi))` times the outer side. This threshold is exact for the universal square-deficit subtraction criterion, not for the geometric inequality. A further conclusion treats an arbitrary parallel rectangular hole whenever both relative side lengths are at least `11/12`, uniformly in their two gap scales and in position. Analytic estimates control all unbounded spectral ranges. Exact integer and rational certificates verify the remaining real parameter intervals, without approximating eigenvalues of a perforated domain.

## Scope

All conclusions concern the ordinary Euclidean Dirichlet Laplacian, with eigenvalues repeated according to multiplicity. The inequality is

`area(Omega) * lambda_k(Omega) >= 4*pi*k`, for every positive integer `k`.

1. **Every parallel square hole:** any size and any strictly interior displacement inside an outer square. This removes the earlier version's centering restriction for intermediate sizes. There is no minimum positive boundary clearance.
2. **Several square holes:** arbitrary positions and orientations, pairwise disjoint and strictly interior, when the sum of side lengths is at most `sqrt(1 - 12/(5*pi))` times the outer side. Sharpness is for the stated scalar subtraction criterion only.
3. **One large parallel rectangular hole:** arbitrary displacement and independent relative side lengths in `[11/12,1)`. The horizontal and vertical gap totals need not be comparable.


## Files

| File | Purpose |
|---|---|
| `r_arun_chandru_Polya_inequalities_displaced_square_holes_perforated_squares.pdf` | Complete manuscript with four appendices |
| `manuscript.tex` | Single self-contained LaTeX source, figure and bibliography included |
| `verify_deficits.py` | Small-hole deficit bounds, partition refinement and analytic constants |
| `verify_corridor.py` | Complete segmented disk census for `8u/9 <= d(u) <= 10u/9`, `u >= 201`, joined to the proved analytic tail |
| `verify_thin_tail.py` | Closed-box scalar certificate and analytic constants for all `q >= 4/5`, `x >= 32` |
| `verify_transport.py` | Independent reconstruction of both finite indexed-eigenvalue certificates |
| `verify_rectangular.py` | Scalar certificates for the large rectangular-hole conclusion |
| `wide_certificate.json`, `thin_certificate.json` | Rational interval and reference-gap **inputs** to the transport verifier |
| `verification_results.json` | Recorded outputs of a complete release replay; not a trusted proof input |
| `SHA256SUMS.txt` | Integrity hashes of the release files |
| `README.md` | This description and reproduction instructions |

The two input JSON files contain explicit rational interval/reference data, not precomputed spectra. The verifier reconstructs every relevant lattice count, cutoff, coefficient, uncovered-index list, and extra gap, and verifies all interval endpoints. No stored PASS flag or certificate-generation program is trusted. The interval references use different denominators (64000 and 192000), which are checked explicitly.

## Reproduce the mathematical certificates

Use Python 3.9 or later. NumPy is the only third-party dependency, needed by the corridor, thin-tail and transport verifiers. The deficit and rectangular verifiers use only Python's standard library. Install NumPy in an environment of your choice if it is not already available:

```text
python -m pip install numpy
```

From this directory run sequentially:

```text
python -B -O verify_deficits.py
python -B -O verify_corridor.py
python -B -O verify_thin_tail.py
python -B -O verify_transport.py
python -B -O verify_rectangular.py
```

Each program prints a JSON summary and exits with a failure if a required check fails. Proof obligations use explicit checks, not removable `assert` statements, and remain active with `-O`. Each program accepts `--output PATH` to save its own result; without this option it does not write result files. Run them sequentially to keep memory use modest. The corridor has a 16 MB frequency segment plus working arrays; other NumPy working arrays also consume memory, so this is not a 16 MB cap on the entire replay process.

The finite comparisons concern full real intervals or boxes. Their infinite-range analytic bounds, variational inequalities, and endpoint arguments are proved in the manuscript, not established by a PASS message. The transport verifier deliberately labels its output as the finite part, whose two high-energy dependencies are verified separately and proved in the paper.

All decisive inequalities use integer or rational arithmetic. Where NumPy floating square roots initialize integer-root proposals, the defining inequalities for those integer roots are checked exactly before use. No floating tolerance determines acceptance. Fixed-width integer ranges are explicitly checked; arbitrary-precision integers handle rational coefficients.

Expected headline results:

- Small-hole deficit census: squared-radius cutoff 1000000; terminal positive-quadrant count 784387; exceptional band 53 handled by the partition argument.
- Strong corridor: squared-radius cutoff 244140625; terminal count 191731890; final failed lower/upper bands 17226/40320, both below `201^2`.
- Thin scalar certificate: 1153748 tested boxes; 865759 accepted leaves; zero unresolved boxes; exact covered integer area 7516192768.
- Transport: 1321 contiguous size intervals; 33926441 interval/index comparisons; 211 additional reference gaps; no uncovered indices.
- Rectangular scalar certificate: 8667 accepted boxes across five comparisons.

The archived output is provided for comparison.

## Compile the manuscript

The source uses `amsart`, standard AMS mathematics packages, Latin Modern, `geometry`, `microtype`, `xurl`, and `hyperref`. Its diagram uses the standard LaTeX picture environment. There are no external figures, bibliography files, included source fragments, shell-escape steps, or local font dependencies beyond those normally supplied by a TeX distribution.

```text
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
```

Alternatively:

```text
latexmk -pdf manuscript.tex
```

or `tectonic manuscript.tex`, which manages its own reruns. LaTeX compilation does not run the mathematical verifiers; run them separately as above.
