# PAD

`pad.exe` evaluates photoelectron angular distributions (PADs) and the
associated anisotropy parameter, beta, from a Dyson-orbital approximation. The
current production path models an outgoing free electron as a linearly
propagating plane wave and is intended for gas-phase photodetachment
simulations where an anion is converted to a neutral molecule plus an electron.

The program reads molecular orbital and basis-set data from a Gaussian FAF
file. In the current workflow, one selected alpha molecular orbital is treated
as the Dyson orbital for a single detachment channel.

## Scientific Model

The default calculation evaluates a length-gauge transition amplitude of the
form

```text
M(k, epsilon) = integral exp(-i k.r) (epsilon.r) psi_D(r) dr
I(theta)      = |M(k, epsilon)|^2
```

where:

- `psi_D(r)` is the selected Dyson orbital, currently supplied as one alpha MO.
- `epsilon` is the electric-field polarization vector.
- `k` is the outgoing photoelectron wave vector.
- `theta` is the angle between `epsilon` and `k`.

For now, global prefactors are not included. The emphasis is on relative
intensities and beta values.

The code keeps the molecular frame fixed and changes the laboratory-frame
vectors. This is the intended path for fixed-orientation calculations and for
future rotational averaging.

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
- Gaussian FAF files containing the needed basis and MO data.

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

## Running

```sh
./pad.exe FAF_FILE MO_INDEX PHOTON_EV BINDING_EV [N_THETA] [N_GRID] [I_PE_TYPE] [LAB_FRAME_TYPE] [N_LAB_THETA] [N_LAB_PHI] [N_CHI]
```

Arguments:

| Argument | Required | Default | Meaning |
| --- | --- | --- | --- |
| `FAF_FILE` | yes | none | Gaussian FAF file to read. |
| `MO_INDEX` | yes | none | Alpha MO index used as the Dyson orbital. |
| `PHOTON_EV` | yes | none | Photon energy in eV. |
| `BINDING_EV` | yes | none | Electron binding/detachment energy in eV. |
| `N_THETA` | no | `5` | Number of theta values from `0` to `pi`. |
| `N_GRID` | no | `101` | Cartesian grid points per axis if FAF grid is absent. |
| `I_PE_TYPE` | no | `0` | Photoelectron model flag. Use `0` for production. |
| `LAB_FRAME_TYPE` | no | `0` | Lab-frame model flag. Use `0` for Cartesian or `1` for sphere-grid. |
| `N_LAB_THETA` | no | `5` | Number of sphere-grid theta points when `LAB_FRAME_TYPE = 1`. |
| `N_LAB_PHI` | no | `8` | Number of sphere-grid phi points when `LAB_FRAME_TYPE = 1`. |
| `N_CHI` | no | `36` | Number of uniform chi points from `0` to `2*pi` used to rotate the PAD scan plane about each `epsilon`. The default gives `10` degree steps. |

The program computes the photoelectron kinetic energy and plane-wave magnitude
from the photon and binding energies:

```text
E_electron(eV)      = PHOTON_EV - BINDING_EV
E_electron(Hartree) = E_electron(eV) / 27.211399
k                   = sqrt(2 * E_electron(Hartree))
```

`PHOTON_EV` must be greater than `BINDING_EV`.

## Examples

Run a plane-wave PAD calculation using an FAF with a stored quadrature grid:

```sh
./pad.exe GTests/006.faf 1 1.100000 1.000000 7 101 0
```

Run a small debug calculation using a fallback Cartesian grid if no FAF grid is
available:

```sh
./pad.exe GTests/001.faf 1 1.100000 1.000000 5 21 0
```

Run the experimental partial-wave path:

```sh
./pad.exe GTests/001.faf 1 1.100000 1.000000 5 21 2
```

Run the sphere-grid lab-frame model with a small orientation grid:

```sh
./pad.exe GTests/006.faf 1 1.100000 1.000000 5 101 0 1 5 8
```

Run the sphere-grid lab-frame model with chi averaging about each polarization
direction:

```sh
./pad.exe GTests/006.faf 1 1.100000 1.000000 5 101 0 1 5 8 8
```

## Output

For each lab-frame orientation, `pad.exe` prints:

- The photoelectron model flag.
- The photon, binding, and resulting kinetic energies.
- The plane-wave magnitude, `k`, in atomic units.
- The lab-frame model flag, number of lab-frame orientations, and orientation
  weight sum.
- The chi quadrature size and chi-weight sum.
- The electric-field polarization vector, `epsilon`.
- The reference perpendicular vector used to define the chi-rotated scan
  planes.
- The theta grid and `I(theta)` values.
- `I(0)`, `I(90)`, and beta from the parallel/perpendicular ratio.
- A fitted beta value from the PAD shape.

The summary table reports the integrated intensity and beta values for each
orientation. When `N_CHI > 1`, the per-orientation PAD and beta values are
first averaged over chi for that `epsilon`. The summary also reports weighted
mean beta values and beta values computed from the weighted
orientation-averaged PAD intensity curve.

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
vectors. The command-line program `pad.f03` is now a thin wrapper around this
driver.

The `pad_options%labFrameType` field controls how the lab-frame vectors are
generated:

```fortran
PAD_LAB_FRAMES_CARTESIAN  ! default 3-axis lab-frame set
PAD_LAB_FRAMES_SPHERE     ! sphere-grid lab-frame set
PAD_LAB_FRAMES_CUSTOM     ! user-supplied vector arrays
```

For the sphere-grid option, set `nLabFrameTheta` and `nLabFramePhi`. For custom
lab frames, fill `labEpsilonVector(3,n)`, `labKPlaneVector(3,n)`, and optionally
`labFrameWeights(n)` and `labFrameLabels(n)`. The lab-frame generator supplies
polarization directions, reference transverse vectors, and outer weights. The
driver then performs a uniform periodic chi average about each `epsilon` using
`nChi` points before combining the weighted `epsilon`-grid results.

The sphere-grid model supplies surface-area weights for the sampled
polarization directions. The chi quadrature supplies equal weights over
`0 <= chi < 2*pi`.

## Interpreting Beta

Two beta estimates are currently printed:

- `beta(ratio)`: computed from `I(0)` and `I(90)`.
- `beta(fit)`: obtained from a least-squares fit to the standard
  `P2(cos(theta))` angular form.

For clean one-channel test cases, these should normally agree closely. Large
disagreements are a useful signal that the angular grid, quadrature, orbital,
or model assumptions need closer inspection.

## Development Notes

Near-term development goals include:

- More robust regression tests for known PAD and beta limits.
- Cleaner handling of FAF metadata and missing quadrature grids.
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

This is currently a lightweight helper test, not a complete validation suite
for the PAD calculation.
