# PAD

`pad.exe` evaluates photoelectron angular distributions (PADs) and the
associated anisotropy parameter, beta, from a Dyson-orbital approximation. The
current production path models an outgoing free electron as a linearly
propagating plane wave and is intended for gas-phase photodetachment
simulations where an anion is converted to a neutral molecule plus an electron.

The program reads coefficient and basis-set data from a Gaussian FAF file. By
default, it expects a DDNO output FAF and uses the labeled DDNO integer metadata
to select `DDNO: ALPHA HOLES` or `DDNO: BETA HOLES` as the Dyson orbital. The
metadata must describe exactly one hole and no particles. Alternatively, the
user can supply `-dyson-mo N`: positive `N` selects column `N` from `ALPHA MO
COEFFICIENTS`, while negative `N` selects column `abs(N)` from `BETA MO
COEFFICIENTS`.

The explicit MO-slot path remains useful when a separate code forms a Dyson
orbital and stores its coefficients in an FAF MO slot for interchange and
visualization. That storage convention does not imply that the selected vector
is an ordinary canonical MO.

## Scientific Model

The default calculation evaluates a length-gauge transition amplitude of the
form

```text
M(k, epsilon) = integral exp(-i k.r) (epsilon.r) psi_D(r) dr
I(theta)      = |M(k, epsilon)|^2
```

where:

- `psi_D(r)` is the selected Dyson orbital, supplied through a DDNO hole record
  or an explicitly selected alpha/beta MO coefficient column.
- `epsilon` is the electric-field polarization vector.
- `k` is the outgoing photoelectron wave vector.
- `theta` is the angle between `epsilon` and `k`.

For now, global prefactors are not included. The emphasis is on relative
intensities and beta values.

The code keeps the molecular frame fixed and changes the laboratory-frame
vectors. This is the intended path for fixed-orientation calculations and
driver-layer rotational averaging.

## Current Status

Use photoelectron model flag `0` for production plane-wave calculations.

Other code paths are present but should be treated as developmental:

- `iPEType = 1`: plane-wave evaluation through the shared
  `dysonMatrixElement` interface.
- `iPEType = 2`: experimental free partial-wave/spherical-Bessel path.

The partial-wave path is useful for development and future decomposition
analysis, but it is not yet validated as a production model.

## Requirements

To build `pad.exe`, you need:

- A Fortran compiler supported by the makefile, currently `gfortran` or
  `nvfortran`.
- MQCPack installed and available through the `mqcinstall` environment
  variable.
- BLAS and LAPACK libraries.
- Gaussian FAF files containing the needed basis data and either compatible
  DDNO records or the explicitly requested alpha/beta MO coefficient array.

For best numerical integration performance, generate FAF files with a stored
Gaussian quadrature grid, for example with Gaussian input options similar to:

```text
output=faf
int(preComputeXCGridPoints)
```

If the FAF does not contain `3D Quadrature Grid`, `pad.exe` falls back to a
simple Cartesian trapezoid grid. That fallback is useful for testing and small
debug cases, but production calculations should prefer a suitable stored
Gaussian quadrature grid.

## Building

From the repository root:

```sh
make pad.exe
```

The makefile expects `mqcinstall` to point to the MQCPack installation. For
example:

```sh
export mqcinstall=/path/to/mqcPack
make pad.exe
```

`make all` builds `pad.exe` and the current lightweight unit-test executables:
`unitTest1.exe`, `unitTest2.exe`, `unitTest3.exe`, and `unitTest4.exe`.

## Test Targets

Run the lightweight unit-style checks with:

```sh
make test-unit
```

This runs `unitTest1.exe`, several `unitTest2.exe` lab-frame setup checks,
`unitTest3.exe` for all 19 Lebedev-Laikov rules, point-count selection, and
spherical product quadrature utilities, and
`unitTest4.exe` for signed-MO and DDNO Dyson-orbital source selection.

Run both the unit-style checks and the PAD regression set with:

```sh
make test
```

## PAD Regression Tests

A small `pad.exe` regression set uses `GTests/006.faf` and compact reference
summaries in `GTests/PAD_Outputs/`. Run it from the repository root with:

```sh
make test-pad
```

The runner checks that `GTests/006.faf` exists before running. If it is absent,
go to `GTests/` and run `006.gjf` with Gaussian to regenerate `006.faf`.

The current cases cover:

- Cartesian lab-frame sampling with `nChi = 1`
- Cartesian lab-frame sampling with `nChi = 36`
- Small sphere-grid lab-frame sampling with `nChi = 1`
- Small sphere-grid lab-frame sampling with `nChi = 4`
- Small axisymmetric lab-frame sampling with `nChi = 1`
- Small axisymmetric lab-frame sampling with `nChi = 4`
- All four PAD print levels and invalid print-level rejection

The test harness compares stable scientific summary quantities with tolerances
rather than diffing full output files, since full output includes wall times.
After an intentional scientific or output-format change, regenerate references
with:

```sh
make update-pad-refs
```

## Running

```sh
./pad.exe -faf FAF_FILE [-dyson-mo MO_INDEX] -photon-ev PHOTON_EV \
  -binding-ev BINDING_EV [-n-theta N_THETA] [-n-grid N_GRID] \
  [-quad QUADRATURE_TYPE] [-lebedev-points N_ANGULAR] \
  [-pe-type I_PE_TYPE] [-lab-frame LAB_FRAME_TYPE] \
  [-lab-theta N_LAB_THETA] [-lab-phi N_LAB_PHI] \
  [-lab-alignment A] [-n-chi N_CHI] [-lmax L_MAX] [-threads N_THREADS] \
  [-print-level terse|normal|extra|debug]
```

Options may be written as `-option value`, `--option value`, or
`--option=value`.

When `-dyson-mo` is omitted, `pad.exe` reads
`DDNO: METADATA INTEGERS LABELS` and `DDNO: METADATA INTEGERS`, locates the
alpha/beta particle and hole counts by their semantic labels, and requires
exactly one hole with no particles. It then reads either `DDNO: ALPHA HOLES` or
`DDNO: BETA HOLES`. Missing, malformed, multielectron, or attachment metadata
terminates the calculation through `mqc_error`.

Core options:

| Option | Required | Default | Meaning |
| --- | --- | --- | --- |
| `-faf` | yes | none | Gaussian FAF file to read. |
| `-dyson-mo` | no | DDNO hole data | Positive `N` selects alpha MO column `N`; negative `N` selects beta MO column `abs(N)`. Zero is invalid when explicitly supplied. |
| `-photon-ev` | yes | none | Photon energy in eV. |
| `-binding-ev` | yes | none | Electron binding/detachment energy in eV. |

Numerical and model options:

| Option | Default | Accepted values | Meaning |
| --- | --- | --- | --- |
| `-n-theta` | `5` | integer, at least `2` | Number of theta values from `0` to `pi`. |
| `-n-grid` | `101` | integer, at least `2` | Cartesian grid points per axis if FAF grid is absent. |
| `-quad` | `auto` | `faf`, `cartesian`, `cart`, `lebedev`, `leb`, `auto` | Matrix-element quadrature source. `auto` preserves the historical behavior: use the FAF grid when present, otherwise Cartesian fallback. |
| `-lebedev-points` | `26` | positive integer, at most `974` | Requested angular point count for the one-center Lebedev product grid. An unavailable count selects the next larger implemented rule. Requests above `974` are errors. Implemented polynomial exactness degree:point-count pairs are `3:6`, `5:14`, `7:26`, `9:38`, `11:50`, `13:74`, `15:86`, `17:110`, `19:146`, `21:170`, `23:194`, `25:230`, `27:266`, `29:302`, `31:350`, `35:434`, `41:590`, `47:770`, and `53:974`. |
| `-pe-type` | `0` | `0`, `1`, `2` | Photoelectron model flag. Use `0` for production. |
| `-lab-frame` | `cartesian` | `cartesian`, `sphere`, `axisymmetric`, `0`, `1`, `2` | Built-in lab-frame model. |
| `-lab-theta` | `5` | integer, at least `3` for `sphere` or `axisymmetric` | Number of sphere-grid theta points when `-lab-frame sphere` or `axisymmetric`. |
| `-lab-phi` | `8` | positive integer | Number of sphere-grid phi points when `-lab-frame sphere` or `axisymmetric`. |
| `-lab-alignment` | `0.0` | real from `-1` to `2` | Axisymmetric alignment strength around lab `z`; internally this is the coefficient `A` in `1 + A P2(cos(theta))`. |
| `-n-chi` | `36` | positive integer | Number of uniform chi points from `0` to `2*pi` used to rotate the PAD scan plane about each `epsilon`. The default gives `10` degree steps. |
| `-lmax` | `6` | nonnegative integer | Maximum angular momentum for developmental partial-wave diagnostics. |
| `-threads` | `1` | positive integer | Number of OpenMP threads requested by the CLI wrapper. |
| `-print-level` | `normal` | `terse`, `normal`, `extra`, `debug` | Controls PAD output detail. `terse` prints the final scientific summary; `normal` adds setup, provenance, and timing output; `extra` adds theta tables; `debug` also enables applicable development diagnostics. |

Common aliases:

| Canonical option | Also accepted |
| --- | --- |
| `-faf` | `-faf-file`, `-file` |
| `-dyson-mo` | `-mo`, `-mo-index`, `-dyson-mo-index` |
| `-photon-ev` | `-photon-energy`, `-photon-energy-ev` |
| `-binding-ev` | `-binding-energy`, `-binding-energy-ev` |
| `-n-theta` | `-theta`, `-n-grid-points-theta` |
| `-n-grid` | `-m-grid`, `-n-grid-points-m` |
| `-quad` | `-quadrature`, `-m-quadrature`, `-spatial-quadrature` |
| `-lebedev-points` | `-lebedev-order`, `-angular-points`, `-n-lebedev` |
| `-pe-type` | `-ipe-type`, `-photoelectron-model` |
| `-lab-frame` | `-lab-frame-type` |
| `-lab-theta` | `-n-lab-theta`, `-n-lab-frame-theta` |
| `-lab-phi` | `-n-lab-phi`, `-n-lab-frame-phi` |
| `-lab-alignment` | `-alignment`, `-lab-align`, `-lab-align-p2`, `-lab-alignment-p2`, `-alignment-p2` |
| `-n-chi` | `-chi` |
| `-threads` | `-omp`, `-n-omp` |
| `-print-level` | `-verbosity` |

The older positional form is still accepted for existing scripts:

```sh
./pad.exe FAF_FILE MO_INDEX PHOTON_EV BINDING_EV [N_THETA] [N_GRID] [I_PE_TYPE] [LAB_FRAME_TYPE] [N_LAB_THETA] [N_LAB_PHI] [N_CHI]
```

The positional `MO_INDEX` remains required and must be nonzero. Its sign uses
the same alpha-positive/beta-negative convention as `-dyson-mo`.

For CLI use, the supported built-in lab-frame models are:

| Model | Numeric value | Meaning |
| --- | ---: | --- |
| `cartesian` | `0` | Three Cartesian polarization directions with equal weights. |
| `sphere` | `1` | Sphere-grid polarization directions with surface-area weights. |
| `axisymmetric` | `2` | Sphere-grid polarization directions reweighted by `1 + A P2(cos(theta))` about lab `z` and renormalized to `4*pi`. |

The programmatic custom model is `PAD_LAB_FRAMES_CUSTOM = -1`. It is intended
for callers that fill `pad_options%labEpsilonVector`,
`pad_options%labKPlaneVector`, and optional custom weights/labels before
calling `runPADCalculation(...)`.

The program computes the photoelectron kinetic energy and plane-wave magnitude
from the photon and binding energies:

```text
E_electron(eV)      = PHOTON_EV - BINDING_EV
E_electron(Hartree) = E_electron(eV) / 27.211399
k                   = sqrt(2 * E_electron(Hartree))
```

`PHOTON_EV` must be greater than `BINDING_EV`.

## Examples

Run from a compatible one-hole DDNO output FAF, allowing `pad.exe` to select
the alpha or beta hole record from the labeled metadata:

```sh
./pad.exe -faf GTests/cn-.ddno.faf -photon-ev 1.100000 \
  -binding-ev 1.000000 -n-theta 5 -n-grid 21 -quad cart -pe-type 0
```

Select beta MO column 2 explicitly:

```sh
./pad.exe -faf FAF_FILE -dyson-mo -2 -photon-ev 1.100000 \
  -binding-ev 1.000000 -n-theta 5 -n-grid 21 -quad cart -pe-type 0
```

Run a plane-wave PAD calculation using an FAF with a stored quadrature grid:

```sh
./pad.exe -faf GTests/006.faf -dyson-mo 1 -photon-ev 1.100000 \
  -binding-ev 1.000000 -n-theta 7 -n-grid 101 -quad faf -pe-type 0
```

Run a small debug calculation using the Cartesian trapezoid grid:

```sh
./pad.exe -faf GTests/001.faf -dyson-mo 1 -photon-ev 1.100000 \
  -binding-ev 1.000000 -n-theta 5 -n-grid 21 -quad cart -pe-type 0
```

Run a calculation using the one-center Lebedev product grid. The angular
point count is independently selectable and is combined with the current
Gauss-Legendre radial default:

```sh
./pad.exe -faf GTests/006.faf -dyson-mo 1 -photon-ev 1.100000 \
  -binding-ev 1.000000 -n-theta 5 -quad leb \
  -lebedev-points 302 -pe-type 0
```

The angular and radial components are separate. The product-grid API also
supports a trapezoid radial rule for tests, and additional radial rules can be
added without changing the Lebedev tables or matrix-element kernels.
The implemented angular point counts are `6`, `14`, `26`, `38`, `50`, `74`,
`86`, `110`, `146`, `170`, `194`, `230`, `266`, `302`, `350`, `434`, `590`,
`770`, and `974`. For example, requesting `300` points selects the 302-point
rule.

Programmatic callers can use
`lebedevSelectPointCount(requestedPoints,roundUp,capAtMaximum)`. The optional
`roundUp` argument defaults to `.true.`; setting it to `.false.` selects the
next smaller rule. An above-maximum request is an error unless the optional
`capAtMaximum` argument is `.true.`, which selects the 974-point rule.

Run the experimental partial-wave path:

```sh
./pad.exe -faf GTests/001.faf -dyson-mo 1 -photon-ev 1.100000 \
  -binding-ev 1.000000 -n-theta 5 -n-grid 21 -pe-type 2
```

Run the sphere-grid lab-frame model with a small orientation grid:

```sh
./pad.exe -faf GTests/006.faf -dyson-mo 1 -photon-ev 1.100000 \
  -binding-ev 1.000000 -n-theta 5 -n-grid 101 -pe-type 0 \
  -lab-frame sphere -lab-theta 5 -lab-phi 8
```

Run the sphere-grid lab-frame model with chi averaging about each polarization
direction:

```sh
./pad.exe -faf GTests/006.faf -dyson-mo 1 -photon-ev 1.100000 \
  -binding-ev 1.000000 -n-theta 5 -n-grid 101 -pe-type 0 \
  -lab-frame sphere -lab-theta 5 -lab-phi 8 -n-chi 8
```

Run the axisymmetric lab-frame model with positive alignment along lab `z`:

```sh
./pad.exe -faf GTests/006.faf -dyson-mo 1 -photon-ev 1.100000 \
  -binding-ev 1.000000 -n-theta 5 -n-grid 101 -pe-type 0 \
  -lab-frame axisymmetric -lab-theta 5 -lab-phi 8 -lab-alignment 0.5
```

## Output

At the default `normal` print level, `pad.exe` reports:

- The photoelectron model flag.
- The photon, binding, and resulting kinetic energies.
- The plane-wave magnitude, `k`, in atomic units.
- The lab-frame model flag, number of lab-frame orientations, and orientation
  weight sum.
- The chi quadrature size and chi-weight sum.
- The selected Dyson-orbital FAF record and MO column, when applicable.
- The four DDNO alpha/beta particle/hole counts for DDNO-driven calculations.
- The electric-field polarization vector, `epsilon`.
- The reference perpendicular vector used to define the chi-rotated scan
  planes.
- Setup and matrix-element timing information.
- The final per-orientation and orientation-averaged scientific summary.

The `terse` level prints only the final scientific summary. The `extra` level
adds the theta grid, `I(theta)` values, `I(0)`, `I(90)`, and beta from the
parallel/perpendicular ratio for each orientation. The `debug` level inherits
all lower-level output and enables applicable development diagnostics.

The summary table reports the integrated intensity and beta values for each
orientation. When `N_CHI > 1`, the per-orientation PAD and beta values are
first averaged over chi for that `epsilon`. The summary also reports weighted
mean beta values and beta values computed from the weighted
orientation-averaged PAD intensity curve.

The inner chi average and outer lab-frame average perform different angular
operations. Chi samples the outgoing-electron azimuth around one fixed
polarization direction. The outer lab-frame weights sample polarization
directions relative to the fixed molecular frame.

Two different integrated-intensity summaries are now printed:

- `theta-integrated intensity`: the PAD integrated only over `theta`
- `solid-angle integrated intensity`: the PAD integrated with the full angular
  measure `sin(theta) dtheta dchi`

## Programmatic Use

The reusable PAD driver lives in `pad_mod.f03`. Other Fortran codes can load or
construct a Gaussian FAF object, fill a `pad_options` object, and call:

```fortran
call runPADCalculation(faf,options,results)
```

where `results` is a `pad_results` object containing the computed kinetic
energy, `k`, theta and chi grids, PAD intensities, theta-integrated
intensities, beta values, solid-angle integrated intensities, and lab-frame
vectors. It also records `dysonSourceLabel`, `dysonFromDDNO`, and the four
`ddnoParticleHoleCounts` in alpha-particle, beta-particle, alpha-hole,
beta-hole order. The command-line program `pad.f03` is now a thin wrapper
around this driver.

For programmatic source selection, `pad_options%dysonMOIndex=0` requests DDNO
metadata (the default), a positive value selects an alpha MO column, and a
negative value selects the corresponding beta MO column.

The `pad_options%labFrameType` field controls how the lab-frame vectors are
generated:

```fortran
PAD_LAB_FRAMES_CUSTOM     ! -1, user-supplied vector arrays
PAD_LAB_FRAMES_CARTESIAN  !  0, default 3-axis lab-frame set
PAD_LAB_FRAMES_SPHERE     !  1, sphere-grid lab-frame set
PAD_LAB_FRAMES_AXISYMMETRIC ! 2, weighted sphere-grid lab-frame set
```

For the sphere-grid and axisymmetric options, set `nLabFrameTheta` and
`nLabFramePhi`. The axisymmetric option also uses `labFrameAlignment`, the
alignment strength around lab `z`. In the current model this is the coefficient
`A` in `1 + A P2(cos(theta))`; valid values are `-1 <= A <= 2`. For custom lab
frames, fill `labEpsilonVector(3,n)`,
`labKPlaneVector(3,n)`, and optionally `labFrameWeights(n)` and
`labFrameLabels(n)`. The lab-frame generator supplies polarization directions,
reference transverse vectors, and outer weights. The driver then performs a
uniform periodic chi average about each `epsilon` using `nChi` points before
combining the weighted `epsilon`-grid results.

The sphere-grid model supplies surface-area weights for the sampled
polarization directions. The axisymmetric model reweights those same
directions and renormalizes the weights to `4*pi`. The chi quadrature supplies
equal weights over `0 <= chi < 2*pi`.

## Interpreting Beta

Two beta estimates are currently printed:

- `beta(ratio)`: computed from `I(0)` and `I(90)`.
- `beta(fit)`: obtained from a least-squares fit to the standard
  `P2(cos(theta))` angular form.

For orientation-averaged results, `averageBetaParaPerp` and `averageBetaFit`
are evaluated from the weighted averaged PAD and are the primary quantities for
comparison with a beta extracted from an averaged experimental PAD.
`meanBetaParaPerp` and `meanBetaFit` instead average beta values that were first
computed separately for each orientation. Those mean values are useful
diagnostics, but they are not generally equal to beta from the averaged PAD
because beta is obtained from ratios of intensities.

The recommended primary reported value is `averageBetaFit`, printed as
`Beta from averaged PAD(fit)`, provided its reported `R**2` is close to one and
the theta grid is converged. `Beta from averaged PAD(ratio)` is an independent
two-angle consistency check. For clean one-channel test cases, these should
normally agree closely. A poor `R**2` or a large disagreement is a signal that
the theta grid, quadrature, orbital, numerical noise, or the assumed
`1 + beta*P2(cos(theta))` angular form needs closer inspection.

## Development Notes

Near-term development goals include:

- More robust regression tests for known PAD and beta limits.
- Cleaner handling of FAF metadata and missing quadrature grids.
- Comparing convergence of the selectable one-center Lebedev product grids
  against Gaussian FAF and Cartesian grids.
- Rotational averaging with a fixed molecular frame and rotated lab frame.
- Validation and extension of the free partial-wave path.
- Later continuum models suitable for photoionization of neutrals, where the
  outgoing electron sees a charged ion rather than a neutral target.
- Performance improvements, especially around quadrature loops and parallelism.

## Related Tests

`unitTest1.exe` exercises spherical-harmonic normalization and spherical Bessel
helper routines:

```sh
make unitTest1.exe
./unitTest1.exe
```

`unitTest2.exe` exercises PAD lab-frame setup behavior:

```sh
make unitTest2.exe
./unitTest2.exe
./unitTest2.exe 1 5 8
./unitTest2.exe 2 5 8 0.5
./unitTest2.exe -1
```

The default run checks the Cartesian model, `1 5 8` checks a small sphere-grid
model, `2 5 8 0.5` checks a positively aligned axisymmetric model, and `-1`
checks the programmatic custom lab-frame path. These are lightweight helper
tests, not a complete validation suite for the PAD calculation.

`unitTest3.exe` exercises every Lebedev-Laikov angular rule from 6 points
(degree 3) through 974 points (degree 53), including analytical Cartesian-
monomial exactness checks. It also tests the one-center product grid with both
trapezoid and Gauss-Legendre radial rules plus helper formulas for choosing
`rMaxBohr` from the AO primitive tail and `nRadial` from the plane-wave
wavelength:

```sh
make unitTest3.exe
./unitTest3.exe
```

`unitTest4.exe` exercises DDNO metadata-label lookup, one-hole record
selection, signed alpha/beta MO loading, and expected rejection of
multielectron or attachment metadata. If `GTests/cn-.ddno.faf` is present, it
also checks the generated DDNO FAF directly:

```sh
make unitTest4.exe
./unitTest4.exe
```
