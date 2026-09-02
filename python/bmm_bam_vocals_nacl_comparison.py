#!/usr/bin/env python3
"""Run BMM parcel-model benchmarks and add BMM to the existing BAM sweep plots.

Intended repository layout
--------------------------
<repo>/
    bam/
        main.exe
        python/
            bam_vocals_nacl_sweeps.py
    bmm/
        main.exe
        namelist.in
        python/
            bmm_bam_vocals_nacl_comparison.py

Run from the repository root, e.g.

    python3 bmm/python/bmm_bam_vocals_nacl_comparison.py

The script imports the BAM sweep module so that VOCALS aerosol, NaCl spray PSD,
thermodynamic conditions, salt-spray choice, and BAM output filenames remain in
sync with the BAM comparison.

BMM experiments
---------------
1. VOCALS aerosol:
   50 logarithmically spaced updrafts from 0.01 to 10 m s-1.

2. VOCALS + NaCl spray:
   the same NaCl mass-mixing-ratio sweep used by BAM, at the BAM salt-sweep
   updraft (normally 0.3 m s-1).

Every BMM parcel starts at RH=0.95, with pinit and tinit taken from the BAM
script, and ascends 500 m.  The comparison is deliberately warm-cloud,
collision-free and uses the full-moving bin scheme:

    ice_flag       = 0
    bin_scheme_flag = 0
    sce_flag        = 0
    kappa_flag      = 0
    sv_flag         = 0
    updraft_type    = 1
    adiabatic_prof  = .true.

The final finite value of NetCDF variable ``ndrop`` is used as the BMM result.

BMM/BAM number units
--------------------
Both BMM and BAM use aerosol number mixing ratio in particles kg-1 dry air and
report activated cloud-drop number in drops kg-1 dry air.  Aerosol numbers are
therefore passed directly from the BAM definitions into BMM with no air-density
conversion.

The script still checks the NetCDF ``ndrop`` units attribute, when present, and
stops if it explicitly reports a volumetric unit.
"""

from __future__ import annotations

import csv
import getpass
import importlib.util
import math
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
BMM_ROOT = SCRIPT_DIR.parent
REPO_ROOT = BMM_ROOT.parent
BAM_ROOT = REPO_ROOT / "bam"

BMM_EXE = BMM_ROOT / "main.exe"
PREFERRED_BMM_NAMELIST = BMM_ROOT / "namelist.in"
PREFERRED_BAM_SCRIPT = BAM_ROOT / "python" / "bam_vocals_nacl_sweeps.py"

USER = getpass.getuser()
BMM_OUTPUT_DIR = Path("/tmp") / USER / "bmm_sweep_output"

# Run the BAM sweep automatically only when its CSV files are missing.
RUN_BAM_IF_MISSING = True

SHOW_PLOTS = True
KEEP_LAST_NETCDF = False


# -----------------------------------------------------------------------------
# BMM comparison settings
# -----------------------------------------------------------------------------
BMM_ASCENT = 500.0                  # m
BMM_ZINIT = 1000.0                  # m
BMM_RHINIT = 0.95                   # fraction

# The BMM updraft sweep uses fewer points because the bin model is much more
# expensive than BAM.
BMM_UPDRAFTS = np.logspace(-2.0, 1.0, 50)

# Limit both time step and vertical displacement per outer BMM step.
BMM_DT_MAX = 10.0                   # s
BMM_DZ_PER_STEP_MAX = 5.0           # m
BMM_DT_MIN = 0.05                   # s

# Warm-cloud activation benchmark.
BMM_N_BINS = 60
BMM_DMINA = 1.0e-9                  # m; below the smallest MCB median diameter
BMM_DMAXA = 3.0e-6                  # m
BMM_SV_FLAG = 0

# BMM and BAM both use aerosol number mixing ratio in kg-1 dry air.
BMM_NUMBER_INPUT_UNITS = "kg-1"

# BMM chemistry: component 1 = ammonium-sulfate-like VOCALS core;
# component 2 = NaCl. Components 3-4 are unused placeholders.
N_COMPS = 4
AS_MOLW = 132.14e-3                 # kg mol-1
AS_RHO = 1770.0                     # kg m-3
AS_NU = 3.0

NACL_MOLW = 58.44e-3                # kg mol-1
NACL_RHO = 2165.0                   # kg m-3
NACL_NU = 2.0


# -----------------------------------------------------------------------------
# Locate/import companion files
# -----------------------------------------------------------------------------
def locate_bam_script() -> Path:
    if PREFERRED_BAM_SCRIPT.exists():
        return PREFERRED_BAM_SCRIPT

    candidates = sorted(
        (BAM_ROOT / "python").glob("bam_vocals_nacl_sweeps*.py")
    )
    if len(candidates) == 1:
        return candidates[0]
    if candidates:
        names = "\n".join(f"  {p}" for p in candidates)
        raise FileNotFoundError(
            "Could not choose the BAM sweep script unambiguously. "
            "Set PREFERRED_BAM_SCRIPT near the top of this file.\n"
            f"Candidates:\n{names}"
        )
    raise FileNotFoundError(
        f"Cannot find a BAM sweep script under {BAM_ROOT / 'python'}"
    )


def locate_bmm_namelist() -> Path:
    if PREFERRED_BMM_NAMELIST.exists():
        return PREFERRED_BMM_NAMELIST

    candidates = sorted(BMM_ROOT.glob("namelist*.in"))
    if len(candidates) == 1:
        return candidates[0]
    if candidates:
        names = "\n".join(f"  {p}" for p in candidates)
        raise FileNotFoundError(
            "Could not choose the BMM base namelist unambiguously. "
            "Set PREFERRED_BMM_NAMELIST near the top of this file.\n"
            f"Candidates:\n{names}"
        )
    raise FileNotFoundError(f"Cannot find a BMM namelist under {BMM_ROOT}")


def import_bam_sweep(path: Path):
    spec = importlib.util.spec_from_file_location("bam_sweep_config", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not import BAM sweep script: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


# -----------------------------------------------------------------------------
# Thermodynamic/unit helpers
# -----------------------------------------------------------------------------
R_DRY_AIR = 287.05  # J kg-1 K-1


def saturation_vapour_pressure_water(t_k: float) -> float:
    """Buck saturation vapour pressure over liquid water [Pa]."""
    tc = t_k - 273.15
    return 611.21 * math.exp((18.678 - tc / 234.5) * (tc / (257.14 + tc)))


def dry_air_density(p_pa: float, t_k: float, rh: float) -> float:
    """Dry-air mass per total moist-air volume [kg dry air m-3]."""
    e = rh * saturation_vapour_pressure_water(t_k)
    if e >= p_pa:
        raise ValueError("Water vapour pressure is not below total pressure")
    return (p_pa - e) / (R_DRY_AIR * t_k)


def number_perkg_to_bmm_input(n_perkg, p_pa, t_k):
    """Return BAM number mixing ratio unchanged for BMM.

    Both models use particles kg-1 dry air.  p_pa and t_k are retained in the
    function signature because the caller also uses them for the parcel setup.
    """
    del p_pa, t_k
    return np.asarray(n_perkg, dtype=float).copy()


def choose_bmm_dt(w: float) -> float:
    if w <= 0.0:
        raise ValueError("BMM comparison requires w > 0")
    return max(
        BMM_DT_MIN,
        min(BMM_DT_MAX, BMM_DZ_PER_STEP_MAX / w),
    )


# -----------------------------------------------------------------------------
# Namelist editing
# -----------------------------------------------------------------------------
def replace_scalar_line(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"(?im)^(\s*{re.escape(key)}\s*=).*$")
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one scalar assignment for {key!r}; "
            f"found {len(matches)}"
        )
    return pattern.sub(lambda m: f"{m.group(1)} {value},", text, count=1)


def replace_namelist_group(text: str, group_name: str, new_group: str) -> str:
    start_match = re.search(
        rf"(?im)^\s*&{re.escape(group_name)}\b", text
    )
    if start_match is None:
        raise RuntimeError(f"Cannot find namelist group &{group_name}")

    next_match = re.search(
        r"(?im)^\s*&[A-Za-z_][A-Za-z0-9_]*\b",
        text[start_match.end():],
    )
    if next_match is None:
        end = len(text)
    else:
        end = start_match.end() + next_match.start()

    return text[:start_match.start()] + new_group.rstrip() + "\n" + text[end:]


def fmt_values(values) -> str:
    return ", ".join(f"{float(v):.16e}" for v in np.asarray(values).ravel())


def pad_internal_modes(values, n_intern: int, fill: float):
    values = np.asarray(values, dtype=float)
    if len(values) > n_intern:
        raise ValueError(
            f"PSD has {len(values)} submodes but BMM n_intern={n_intern}"
        )
    out = np.full(n_intern, fill, dtype=float)
    out[: len(values)] = values
    return out


def make_bmm_modal_arrays(
    bam,
    p_test: float,
    t_test: float,
    salt_n_perkg=None,
    salt_d=None,
    salt_logsig=None,
):
    """Build BMM n_aer1/d_aer1/sig_aer1 arrays.

    BMM represents the three VOCALS lognormals as internal submodes of one
    external chemical mode.  Added NaCl is a second external chemical mode.
    """
    n_intern = len(np.asarray(bam.VOCALS_N))
    if n_intern != 3:
        raise ValueError(
            "This script expects the current three-submode VOCALS definition"
        )

    ext_n_perkg = [np.asarray(bam.VOCALS_N, dtype=float)]
    ext_d = [np.asarray(bam.VOCALS_D, dtype=float)]
    ext_sig = [np.asarray(bam.VOCALS_LOGSIG, dtype=float)]

    if salt_n_perkg is not None:
        salt_n_perkg = pad_internal_modes(salt_n_perkg, n_intern, 0.0)
        salt_d = pad_internal_modes(salt_d, n_intern, 100.0e-9)
        salt_logsig = pad_internal_modes(salt_logsig, n_intern, 0.3)
        ext_n_perkg.append(salt_n_perkg)
        ext_d.append(salt_d)
        ext_sig.append(salt_logsig)

    n_perkg = np.column_stack(ext_n_perkg)   # (n_intern, n_mode)
    d = np.column_stack(ext_d)
    sig = np.column_stack(ext_sig)

    # BMM and BAM use the same particle-number mixing-ratio basis (kg-1).
    n_input = number_perkg_to_bmm_input(n_perkg, p_test, t_test)
    return n_input, d, sig


def make_aerosol_setup_group(n_mode: int, n_intern: int) -> str:
    return f"""&aerosol_setup
    n_mode        = {n_mode},
    n_intern      = {n_intern},
    n_sv          = 10,
    sv_flag       = {BMM_SV_FLAG},
    n_bins        = {BMM_N_BINS},
    n_comps       = {N_COMPS},
    n_inp_classes = 0
/"""


def make_aerosol_spec_group(n_input, d, sig) -> str:
    n_input = np.asarray(n_input, dtype=float)
    d = np.asarray(d, dtype=float)
    sig = np.asarray(sig, dtype=float)

    if not (n_input.shape == d.shape == sig.shape):
        raise ValueError("BMM modal number/diameter/width arrays must match")

    n_intern, n_mode = n_input.shape

    # mass_frac_aer1 is indexed (external mode, component).
    mass_frac = np.zeros((n_mode, N_COMPS), dtype=float)
    mass_frac[0, 0] = 1.0                 # VOCALS ammonium sulfate
    if n_mode > 1:
        mass_frac[1, 1] = 1.0             # injected NaCl

    # Unused components remain chemically valid placeholders.
    molw = np.array([AS_MOLW, NACL_MOLW, AS_MOLW, AS_MOLW])
    rho = np.array([AS_RHO, NACL_RHO, AS_RHO, AS_RHO])
    nu = np.array([AS_NU, NACL_NU, AS_NU, AS_NU])

    # kappa is not used when kappa_flag=0, but set it explicitly to zero for
    # the standard molecular Koehler comparison requested here.
    kappa = np.zeros(N_COMPS)

    # Arrays are flattened in Fortran storage order to match the explicit
    # subarray shapes used in the supplied namelist.
    n_flat = n_input.flatten(order="F")
    d_flat = d.flatten(order="F")
    sig_flat = sig.flatten(order="F")
    mf_flat = mass_frac.flatten(order="F")

    return f"""&aerosol_spec
    n_aer1(1:{n_intern},1:{n_mode}) = {fmt_values(n_flat)},
    d_aer1(1:{n_intern},1:{n_mode}) = {fmt_values(d_flat)},
    sig_aer1(1:{n_intern},1:{n_mode}) = {fmt_values(sig_flat)},

    dmina = {BMM_DMINA:.16e},
    dmaxa = {BMM_DMAXA:.16e},

    mass_frac_aer1(1:{n_mode},1:{N_COMPS}) = {fmt_values(mf_flat)},

    molw_core1(1:{N_COMPS}) = {fmt_values(molw)},
    density_core1(1:{N_COMPS}) = {fmt_values(rho)},
    nu_core1(1:{N_COMPS}) = {fmt_values(nu)},
    kappa_core1(1:{N_COMPS}) = {fmt_values(kappa)},

    inp_category(1:{N_COMPS}) = 'none', 'none', 'none', 'none',

    org_content1(1:10) = 0.005, 0.01, 0.02, 0.03, 0.06, 0.08, 0.16, 0.3, 0.42, 0.8,
    molw_org1(1:10) = 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
    kappa_org1(1:10) = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    density_org1(1:10) = 1500., 1500., 1500., 1500., 1500., 1500., 1500., 1500., 1500., 1500.,
    delta_h_vap1(1:10) = 150., 150., 150., 150., 150., 150., 150., 150., 150., 150.,
    nu_org1(1:10) = 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
    log_c_star1(1:10) = -6., -5., -4., -3., -2., -1., 0., 1., 2., 3.
/"""


def write_bmm_namelist(
    base_text: str,
    path: Path,
    output_nc: Path,
    bam,
    w: float,
    salt_n_perkg=None,
    salt_d=None,
    salt_logsig=None,
):
    p_test = float(bam.P_TEST)
    t_test = float(bam.T_TEST)
    dt = choose_bmm_dt(w)
    runtime = BMM_ASCENT / w + 2.0 * dt

    n_input, d, sig = make_bmm_modal_arrays(
        bam,
        p_test,
        t_test,
        salt_n_perkg=salt_n_perkg,
        salt_d=salt_d,
        salt_logsig=salt_logsig,
    )
    n_intern, n_mode = n_input.shape

    text = base_text

    # Thermodynamic/ascent setup.
    replacements = {
        "outputfile": f"'{output_nc}'",
        "runtime": f"{runtime:.16e}",
        "dt": f"{dt:.16e}",
        "zinit": f"{BMM_ZINIT:.16e}",
        "tpert": "0.0",
        "use_prof_for_tprh": ".false.",
        "winit": f"{w:.16e}",
        "tinit": f"{t_test:.16e}",
        "pinit": f"{p_test:.16e}",
        "rhinit": f"{BMM_RHINIT:.16e}",
        "ice_flag": "0",
        "bin_scheme_flag": "0",
        "sce_flag": "0",
        "kappa_flag": "0",
        "updraft_type": "1",
        "adiabatic_prof": ".true.",
        "entrain_period": "0",
        "ent_rate": "0.0",
        "fallout_flag": ".false.",
        "z_ctop": f"{BMM_ZINIT + BMM_ASCENT:.16e}",
    }
    for key, value in replacements.items():
        text = replace_scalar_line(text, key, value)

    text = replace_namelist_group(
        text,
        "aerosol_setup",
        make_aerosol_setup_group(n_mode, n_intern),
    )
    text = replace_namelist_group(
        text,
        "aerosol_spec",
        make_aerosol_spec_group(n_input, d, sig),
    )

    path.write_text(text)


# -----------------------------------------------------------------------------
# NetCDF output
# -----------------------------------------------------------------------------
def _units_are_volumetric(units: str) -> bool:
    u = units.lower().replace(" ", "")
    markers = (
        "m-3",
        "m^-3",
        "m**-3",
        "/m3",
        "m^{-3}",
    )
    return any(marker in u for marker in markers)


def _last_finite_value(array) -> float:
    arr = np.asarray(array, dtype=float).ravel()
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        raise RuntimeError("ndrop contains no finite values")
    return float(finite[-1])


def read_final_ndrop(path: Path) -> tuple[float, str]:
    """Read the final finite ndrop value and its units string if available."""
    try:
        from netCDF4 import Dataset
    except ImportError:
        Dataset = None

    if Dataset is not None:
        with Dataset(path, "r") as ds:
            if "ndrop" not in ds.variables:
                raise KeyError(f"'ndrop' is not present in {path}")
            var = ds.variables["ndrop"]
            units = str(getattr(var, "units", ""))
            data = np.ma.filled(var[:], np.nan)
            value = _last_finite_value(data)
    else:
        try:
            import xarray as xr
        except ImportError as exc:
            raise ImportError(
                "Reading BMM output requires either netCDF4 or xarray"
            ) from exc

        with xr.open_dataset(path) as ds:
            if "ndrop" not in ds:
                raise KeyError(f"'ndrop' is not present in {path}")
            units = str(ds["ndrop"].attrs.get("units", ""))
            value = _last_finite_value(ds["ndrop"].values)

    if units and _units_are_volumetric(units):
        raise RuntimeError(
            f"BMM ndrop reports volumetric units {units!r}, but BAM plots "
            "drop number in kg-1 dry air. Add an explicit BMM output-unit "
            "conversion before comparing these curves."
        )

    return value, units


# -----------------------------------------------------------------------------
# Model execution and resumable sweeps
# -----------------------------------------------------------------------------
def run_bmm(
    base_text: str,
    case_nml: Path,
    output_nc: Path,
    bam,
    w: float,
    salt_n_perkg=None,
    salt_d=None,
    salt_logsig=None,
) -> float:
    if output_nc.exists():
        output_nc.unlink()

    write_bmm_namelist(
        base_text,
        case_nml,
        output_nc,
        bam,
        w,
        salt_n_perkg=salt_n_perkg,
        salt_d=salt_d,
        salt_logsig=salt_logsig,
    )

    proc = subprocess.run(
        [str(BMM_EXE), str(case_nml)],
        cwd=str(BMM_ROOT),
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"BMM failed for w={w:g} m s-1\n"
            f"Generated namelist: {case_nml}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )
    if not output_nc.exists():
        raise FileNotFoundError(
            f"BMM finished without creating expected output {output_nc}"
        )

    ndrop, units = read_final_ndrop(output_nc)
    if units:
        print(f"      ndrop={ndrop:.6e} [{units}]")
    else:
        print(f"      ndrop={ndrop:.6e} [units attribute absent; assumed kg-1]")

    if not KEEP_LAST_NETCDF:
        output_nc.unlink(missing_ok=True)

    return ndrop


def load_or_create_cache(path: Path, x: np.ndarray) -> np.ndarray:
    y = np.full_like(x, np.nan, dtype=float)
    if not path.exists():
        return y

    try:
        old = np.loadtxt(path, delimiter=",", skiprows=1)
        old = np.atleast_2d(old)
        if old.shape == (len(x), 2) and np.allclose(
            old[:, 0], x, rtol=1e-12, atol=0.0
        ):
            y[:] = old[:, 1]
            print(f"  Resuming from {path}")
    except Exception:
        pass
    return y


def save_cache(path: Path, xname: str, x: np.ndarray, y: np.ndarray):
    np.savetxt(
        path,
        np.column_stack((x, y)),
        delimiter=",",
        header=f"{xname},BMM_kg-1",
        comments="",
    )


def bmm_updraft_sweep(base_text, case_nml, output_nc, bam):
    print("Running BMM VOCALS updraft sweep...")
    csv_path = BMM_OUTPUT_DIR / "bmm_vocals_updraft_sweep.csv"
    ndrop = load_or_create_cache(csv_path, BMM_UPDRAFTS)

    for j, w in enumerate(BMM_UPDRAFTS):
        if np.isfinite(ndrop[j]):
            continue
        print(f"  [{j+1:02d}/{len(BMM_UPDRAFTS)}] w={w:.6g} m s-1")
        ndrop[j] = run_bmm(
            base_text,
            case_nml,
            output_nc,
            bam,
            float(w),
        )
        save_cache(csv_path, "w_m_s-1", BMM_UPDRAFTS, ndrop)

    return BMM_UPDRAFTS, ndrop


def bmm_salt_sweep(base_text, case_nml, output_nc, bam):
    weights, salt_logsig, salt_dm, spray_name = bam.get_spray_psd(
        bam.spray_method
    )
    q_values = np.asarray(bam.NACL_MR, dtype=float)
    w = float(bam.SALT_SWEEP_W)
    stem = spray_name.lower()

    print(
        f"Running BMM VOCALS + NaCl sweep using {spray_name} at "
        f"w={w:g} m s-1..."
    )
    csv_path = BMM_OUTPUT_DIR / f"bmm_vocals_nacl_{stem}_sweep.csv"
    ndrop = load_or_create_cache(csv_path, q_values)

    for j, q_nacl in enumerate(q_values):
        if np.isfinite(ndrop[j]):
            continue

        salt_n_perkg = bam.nacl_number_from_mixing_ratio(
            float(q_nacl),
            weights,
            salt_logsig,
            salt_dm,
        )
        print(
            f"  [{j+1:02d}/{len(q_values)}] "
            f"q_NaCl={q_nacl:.6e} kg kg-1"
        )
        ndrop[j] = run_bmm(
            base_text,
            case_nml,
            output_nc,
            bam,
            w,
            salt_n_perkg=salt_n_perkg,
            salt_d=salt_dm,
            salt_logsig=salt_logsig,
        )
        save_cache(csv_path, "NaCl_kg_kg-1", q_values, ndrop)

    return q_values, ndrop, spray_name


# -----------------------------------------------------------------------------
# BAM CSVs and combined plots
# -----------------------------------------------------------------------------
def bam_output_suffix(bam) -> str:
    if hasattr(bam, "output_suffix"):
        return str(bam.output_suffix())
    if hasattr(bam, "sv_suffix"):
        return str(bam.sv_suffix())
    return ""


def expected_bam_csvs(bam):
    suffix = bam_output_suffix(bam)
    weights, salt_logsig, salt_dm, spray_name = bam.get_spray_psd(
        bam.spray_method
    )
    del weights, salt_logsig, salt_dm
    bam_dir = Path(bam.OUTPUT_DIR)
    return (
        bam_dir / f"vocals_updraft_sweep{suffix}.csv",
        bam_dir / f"vocals_nacl_{spray_name.lower()}_sweep{suffix}.csv",
        spray_name,
    )


def ensure_bam_results(bam_script: Path, bam):
    up_csv, salt_csv, _ = expected_bam_csvs(bam)
    if up_csv.exists() and salt_csv.exists():
        return

    if not RUN_BAM_IF_MISSING:
        raise FileNotFoundError(
            "BAM sweep CSV files are missing. Run the BAM sweep first, or set "
            "RUN_BAM_IF_MISSING=True."
        )

    print("BAM sweep CSVs are missing; running the BAM sweep first...")
    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    proc = subprocess.run(
        [sys.executable, str(bam_script)],
        cwd=str(REPO_ROOT),
        text=True,
        capture_output=True,
        env=env,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "BAM sweep failed while generating comparison CSV files.\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )

    if not up_csv.exists() or not salt_csv.exists():
        raise FileNotFoundError(
            "BAM sweep finished, but the expected CSV files were not created:\n"
            f"  {up_csv}\n  {salt_csv}"
        )


def load_bam_table(path: Path):
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    data = np.atleast_2d(data)
    return header, data


def bam_label_from_header(name: str) -> str:
    label = name
    if label.endswith("_kg-1"):
        label = label[:-5]
    return label.replace("_", " ")


def plot_updraft_comparison(bam, bmm_x, bmm_ndrop):
    up_csv, _, _ = expected_bam_csvs(bam)
    header, data = load_bam_table(up_csv)

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for j in range(1, data.shape[1]):
        ax.plot(data[:, 0], data[:, j], label=bam_label_from_header(header[j]))

    ax.plot(
        bmm_x,
        bmm_ndrop,
        marker="o",
        markersize=3,
        linewidth=1.5,
        label="BMM",
    )
    ax.set_xscale("log")
    ax.set_xlabel(r"Updraft speed (m s$^{-1}$)")
    ax.set_ylabel(r"Activated drop number (kg$^{-1}$ dry air)")
    ax.set_title(
        f"VOCALS aerosol activation: BAM vs BMM "
        f"({BMM_ASCENT:g} m BMM ascent)"
    )
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()

    suffix = bam_output_suffix(bam)
    path = BMM_OUTPUT_DIR / f"vocals_updraft_sweep_with_bmm{suffix}.png"
    fig.savefig(path, dpi=200)
    return fig, path


def plot_salt_comparison(bam, bmm_x, bmm_ndrop, spray_name):
    _, salt_csv, _ = expected_bam_csvs(bam)
    header, data = load_bam_table(salt_csv)

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for j in range(1, data.shape[1]):
        ax.plot(data[:, 0], data[:, j], label=bam_label_from_header(header[j]))

    ax.plot(
        bmm_x,
        bmm_ndrop,
        marker="o",
        markersize=3,
        linewidth=1.5,
        label="BMM",
    )
    ax.set_xscale("log")
    ax.set_xlabel(r"NaCl mixing ratio (kg kg$^{-1}$)")
    ax.set_ylabel(r"Activated drop number (kg$^{-1}$ dry air)")
    ax.set_title(
        f"VOCALS + NaCl ({spray_name}): BAM vs BMM, "
        f"w = {float(bam.SALT_SWEEP_W):g} m s$^{{-1}}$"
    )
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.tight_layout()

    suffix = bam_output_suffix(bam)
    path = (
        BMM_OUTPUT_DIR
        / f"vocals_nacl_{spray_name.lower()}_sweep_with_bmm{suffix}.png"
    )
    fig.savefig(path, dpi=200)
    return fig, path


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main():
    if not BMM_EXE.exists():
        raise FileNotFoundError(
            f"Cannot find BMM executable {BMM_EXE}. Build BMM first."
        )

    bam_script = locate_bam_script()
    bmm_namelist = locate_bmm_namelist()
    bam = import_bam_sweep(bam_script)

    required_bam_names = (
        "P_TEST",
        "T_TEST",
        "VOCALS_N",
        "VOCALS_D",
        "VOCALS_LOGSIG",
        "NACL_MR",
        "SALT_SWEEP_W",
        "spray_method",
        "get_spray_psd",
        "nacl_number_from_mixing_ratio",
        "OUTPUT_DIR",
    )
    missing = [name for name in required_bam_names if not hasattr(bam, name)]
    if missing:
        raise AttributeError(
            "The BAM sweep script is missing required definitions: "
            + ", ".join(missing)
        )

    if getattr(bam, "SV_FLAG", 0) != 0 and BMM_SV_FLAG == 0:
        print(
            "WARNING: BAM semi-volatiles are ON but this BMM benchmark has "
            "BMM_SV_FLAG=0, so the curves are not compositionally equivalent."
        )

    BMM_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"BAM sweep module: {bam_script}")
    print(f"BMM base namelist: {bmm_namelist}")
    print(f"BMM executable: {BMM_EXE}")
    print(f"BMM output directory: {BMM_OUTPUT_DIR}")
    print(
        f"BMM initial state: T={float(bam.T_TEST):g} K, "
        f"P={float(bam.P_TEST):g} Pa, RH={BMM_RHINIT:g}"
    )
    print(f"BMM ascent: {BMM_ASCENT:g} m")
    print(
        "BMM aerosol-number input: particles kg-1 dry air "
        "(passed directly from BAM; no density conversion)"
    )

    ensure_bam_results(bam_script, bam)

    base_text = bmm_namelist.read_text()
    case_nml = BMM_OUTPUT_DIR / "_bmm_sweep_case.in"
    output_nc = BMM_OUTPUT_DIR / "_bmm_sweep_case.nc"

    up_x, up_ndrop = bmm_updraft_sweep(
        base_text, case_nml, output_nc, bam
    )
    salt_x, salt_ndrop, spray_name = bmm_salt_sweep(
        base_text, case_nml, output_nc, bam
    )

    fig1, plot1 = plot_updraft_comparison(bam, up_x, up_ndrop)
    fig2, plot2 = plot_salt_comparison(
        bam, salt_x, salt_ndrop, spray_name
    )

    print("\nCombined plots:")
    print(f"  {plot1}")
    print(f"  {plot2}")
    print(f"BMM CSV data are also in {BMM_OUTPUT_DIR}")

    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close(fig1)
        plt.close(fig2)


if __name__ == "__main__":
    main()
