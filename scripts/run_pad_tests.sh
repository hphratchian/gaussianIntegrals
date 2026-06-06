#!/usr/bin/env bash
set -euo pipefail

update_refs=false
if [[ "${1:-}" == "--update" ]]; then
  update_refs=true
elif [[ $# -gt 0 ]]; then
  echo "Usage: $0 [--update]" >&2
  exit 2
fi

faf="GTests/006.faf"
gjf="006.gjf"
ref_dir="GTests/PAD_Outputs"
work_dir="${TMPDIR:-/tmp}/gaussianIntegrals-pad-tests.$$"
mkdir -p "$work_dir"
trap 'rm -rf "$work_dir"' EXIT

if [[ ! -f "$faf" ]]; then
  cat >&2 <<MSG
PAD regression tests require $faf.

Please go to GTests and run $gjf with Gaussian to regenerate 006.faf, then
rerun this test script.
MSG
  exit 2
fi

if [[ ! -x ./pad.exe ]]; then
  echo "PAD regression tests require ./pad.exe. Run 'make pad.exe' first." >&2
  exit 2
fi

cases=(
  "006_cartesian_nchi1|GTests/006.faf 1 1.100000 1.000000 3 11 0 0 5 8 1"
  "006_cartesian_nchi36|-faf GTests/006.faf -dyson-mo 1 -photon-ev 1.100000 -binding-ev 1.000000 -n-theta 3 -n-grid 11 -pe-type 0 -lab-frame cartesian -n-chi 36"
  "006_sphere_small_nchi1|-faf GTests/006.faf -dyson-mo 1 -photon-ev 1.100000 -binding-ev 1.000000 -n-theta 3 -n-grid 11 -pe-type 0 -lab-frame sphere -lab-theta 3 -lab-phi 4 -n-chi 1"
  "006_sphere_small_nchi4|-faf GTests/006.faf -dyson-mo 1 -photon-ev 1.100000 -binding-ev 1.000000 -n-theta 3 -n-grid 11 -pe-type 0 -lab-frame sphere -lab-theta 3 -lab-phi 4 -n-chi 4"
  "006_axisymmetric_small_nchi1|-faf GTests/006.faf -dyson-mo 1 -photon-ev 1.100000 -binding-ev 1.000000 -n-theta 3 -n-grid 11 -pe-type 0 -lab-frame axisymmetric -lab-theta 3 -lab-phi 4 -lab-alignment 0.5 -n-chi 1"
  "006_axisymmetric_small_nchi4|-faf GTests/006.faf -dyson-mo 1 -photon-ev 1.100000 -binding-ev 1.000000 -n-theta 3 -n-grid 11 -pe-type 0 -lab-frame axisymmetric -lab-theta 3 -lab-phi 4 -lab-alignment 0.5 -n-chi 4"
)

extract_summary() {
  awk '
    function trim(x) { gsub(/^[[:space:]]+|[[:space:]]+$/, "", x); return x }
    function first_value_after_equals(line, value, parts) {
      sub(/^[^=]*=[[:space:]]*/, "", line)
      value = line
      sub(/[[:space:]].*$/, "", value)
      return trim(value)
    }
    /Electron kinetic energy=/ {
      line = $0
      sub(/^.*Electron kinetic energy=[[:space:]]*/, "", line)
      split(line, a, "eV")
      photoelectronEnergyEV = trim(a[1])
      line = a[2]
      sub(/^[[:space:]]*=[[:space:]]*/, "", line)
      split(line, b, "Eh")
      photoelectronEnergyHartree = trim(b[1])
    }
    /Photoelectron k[[:space:]]*=/ {
      kMag = first_value_after_equals($0)
    }
    /Lab-frame orientations[[:space:]]*=/ {
      line = $0
      sub(/^.*Lab-frame orientations[[:space:]]*=[[:space:]]*/, "", line)
      split(line, a, "weight sum")
      nLabFrames = trim(a[1])
      line = a[2]
      sub(/^[[:space:]]*=[[:space:]]*/, "", line)
      labFrameWeightSum = trim(line)
    }
    /Chi quadrature points[[:space:]]*=/ {
      line = $0
      sub(/^.*Chi quadrature points[[:space:]]*=[[:space:]]*/, "", line)
      split(line, a, "weight sum")
      nChi = trim(a[1])
      line = a[2]
      sub(/^[[:space:]]*=[[:space:]]*/, "", line)
      chiWeightSum = trim(line)
    }
    /Averaged theta-integrated intensity[[:space:]]*=/ {
      averageThetaIntegratedIntensity = first_value_after_equals($0)
    }
    /Averaged solid-angle integrated intensity[[:space:]]*=/ {
      averageSolidAngleIntegratedIntensity = first_value_after_equals($0)
    }
    /Beta from averaged PAD\(ratio\)[[:space:]]*=/ {
      averageBetaParaPerp = first_value_after_equals($0)
    }
    /Beta from averaged PAD\(fit\)[[:space:]]*=/ {
      line = $0
      sub(/^.*Beta from averaged PAD[(]fit[)][[:space:]]*=[[:space:]]*/, "", line)
      sub(/[[:space:]].*$/, "", line)
      averageBetaFit = trim(line)
    }
    END {
      print "# key value abs_tol rel_tol"
      print "photoelectronEnergyEV", photoelectronEnergyEV, "5.0e-7", "5.0e-7"
      print "photoelectronEnergyHartree", photoelectronEnergyHartree, "5.0e-9", "5.0e-7"
      print "kMag", kMag, "5.0e-7", "5.0e-7"
      print "nLabFrames", nLabFrames, "0.0", "0.0"
      print "labFrameWeightSum", labFrameWeightSum, "5.0e-6", "5.0e-7"
      print "nChi", nChi, "0.0", "0.0"
      print "chiWeightSum", chiWeightSum, "5.0e-6", "5.0e-7"
      print "averageThetaIntegratedIntensity", averageThetaIntegratedIntensity, "5.0e-11", "5.0e-7"
      print "averageSolidAngleIntegratedIntensity", averageSolidAngleIntegratedIntensity, "5.0e-26", "5.0e-7"
      print "averageBetaParaPerp", averageBetaParaPerp, "5.0e-7", "5.0e-7"
      print "averageBetaFit", averageBetaFit, "5.0e-7", "5.0e-7"
    }
  ' "$1"
}

compare_summary() {
  awk '
    function is_number(x) {
      return x ~ /^[-+]?(([0-9]+([.][0-9]*)?)|([.][0-9]+))([Ee][-+]?[0-9]+)?$/
    }
    FNR == NR {
      if($1 !~ /^#/ && NF >= 4) {
        if(!is_number($2) || !is_number($3) || !is_number($4)) {
          printf("  nonnumeric reference value for %s\n", $1)
          failures++
        }
        ref[$1] = $2
        absTol[$1] = $3
        relTol[$1] = $4
      }
      next
    }
    $1 !~ /^#/ && NF >= 4 {
      seen[$1] = 1
      if(!is_number($2) || !is_number($3) || !is_number($4)) {
        printf("  nonnumeric current value for %s\n", $1)
        failures++
        next
      }
      current = $2 + 0.0
      expected = ref[$1] + 0.0
      tolerance = absTol[$1] + relTol[$1] * (expected < 0 ? -expected : expected)
      diff = current - expected
      if(diff < 0) diff = -diff
      if(!($1 in ref)) {
        printf("  unexpected key %s\n", $1)
        failures++
      } else if(diff > tolerance) {
        printf("  %s expected %.12e got %.12e diff %.3e tol %.3e\n",  \
          $1, expected, current, diff, tolerance)
        failures++
      }
    }
    END {
      for(key in ref) {
        if(!(key in seen)) {
          printf("  missing key %s\n", key)
          failures++
        }
      }
      exit failures ? 1 : 0
    }
  ' "$1" "$2"
}

failures=0
for case_spec in "${cases[@]}"; do
  name="${case_spec%%|*}"
  args="${case_spec#*|}"
  out_file="$work_dir/$name.out"
  summary_file="$work_dir/$name.summary"
  ref_file="$ref_dir/$name.ref"

  echo "Running $name"
  # shellcheck disable=SC2086
  ./pad.exe $args > "$out_file"
  extract_summary "$out_file" > "$summary_file"

  if $update_refs; then
    cp "$summary_file" "$ref_file"
    echo "  updated $ref_file"
  else
    if [[ ! -f "$ref_file" ]]; then
      echo "  missing reference $ref_file" >&2
      failures=$((failures+1))
      continue
    fi
    if compare_summary "$ref_file" "$summary_file"; then
      echo "  PASS"
    else
      echo "  FAIL"
      failures=$((failures+1))
    fi
  fi
done

if [[ "$failures" -ne 0 ]]; then
  echo "PAD regression tests failed: $failures case(s)." >&2
  exit 1
fi

if $update_refs; then
  echo "PAD regression references updated."
else
  echo "PAD regression tests passed."
fi
