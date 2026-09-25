#!/usr/bin/env python3
"""Run adiabatic-parcel NaCl sweeps for added-aerosol PSDs 1--6.

The ammonium-sulphate background is fixed to PSD2.  The added aerosol is
NaCl with PSD1--PSD6.  Aerosol number concentrations can optionally be scaled
with updraft speed using

    N_aer / N_aer,ref = (w / w_ref)**(3/2)

as an approximate supersaturation-similarity scaling.

The same model outputs can be analysed in either of two ways:

  * ``ndrop`` : use the model's activated-drop ``ndrop`` field;
  * ``dgt``   : count liquid drops with D > a chosen threshold (2 um by default)
                using ``nwat`` and ``mwat``.

Typical use
-----------
Run all 6 x 50 simulations and plot the default diagnostic::

    python3 allPSDRunsAIDA.py

Replot existing outputs using D > 2 um instead of ndrop::

    python3 allPSDRunsAIDA.py --plot-only --diagnostic dgt --dmin-um 2

The most commonly changed settings are in the Experiment definition section
below, especially ``WINIT`` and ``DROP_DIAGNOSTIC``.
"""

from __future__ import annotations

import argparse
import csv
import getpass
import re
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


# -----------------------------------------------------------------------------
# Experiment definition -- edit these values
# -----------------------------------------------------------------------------

BMM_RUN = False

BACKGROUND_PSD = 2
ADDED_PSDS = (1, 2, 3, 4, 5, 6)
N_RUNS = 50

# Updraft for this experiment [m s^-1].
# Change this to e.g. 0.3 for the lower-updraft similarity experiment.
WINIT = 0.3

# Reference chamber-like experiment used for the aerosol-number scaling.
REFERENCE_UPDRAFT_MS = 1.3
REFERENCE_BACKGROUND_NUMBER_CM3 = 3000.0
REFERENCE_MAX_ADDED_NUMBER_CM3 = 5000.0

# If true, scale BOTH the background and added-aerosol number concentrations as
# (WINIT / REFERENCE_UPDRAFT_MS)**1.5.  Set false to use the reference aerosol
# concentrations regardless of updraft.
SCALE_AEROSOL_WITH_UPDRAFT = True
AEROSOL_SCALING_EXPONENT = 1.0 #1.5

# Drop-number diagnostic used when plotting/analyzing.
#   "ndrop" : model activated-drop diagnostic
#   "dgt"   : N for drops larger than DROP_DIAMETER_THRESHOLD_UM
DROP_DIAGNOSTIC = "ndrop"
DROP_DIAMETER_THRESHOLD_UM = 2.0

# Initial state used to convert prescribed number per m3 to model number per kg.
RHO_AIR_INITIAL = 100000.0 / 280.0 / 287.0
RHO_WATER = 1000.0
R_D = 287.0

KAPPA_BACK = 0.61       # ammonium sulphate
KAPPA_ADD = 1.28        # NaCl
DENSITY_BACK = 1770.0   # kg m^-3
DENSITY_ADD = 2165.0    # kg m^-3


# Derived experiment concentrations.
if SCALE_AEROSOL_WITH_UPDRAFT:
    AEROSOL_NUMBER_SCALE = (
        WINIT / REFERENCE_UPDRAFT_MS
    ) ** AEROSOL_SCALING_EXPONENT
else:
    AEROSOL_NUMBER_SCALE = 1.0

BACKGROUND_NUMBER_CM3 = (
    REFERENCE_BACKGROUND_NUMBER_CM3 * AEROSOL_NUMBER_SCALE
)
MAX_ADDED_NUMBER_CM3 = (
    REFERENCE_MAX_ADDED_NUMBER_CM3 * AEROSOL_NUMBER_SCALE
)


# PSD definitions copied from runsDefineAIDA.py.
# Nfrac gives the relative number in each lognormal component.
PSD = {
    1: {
        "Nfrac": np.array([0.49, 0.38, 1.0e-8]),
        "logSigma": np.array([0.25, 0.84, 0.25]),
        "Dm": np.array([0.247e-6, 0.205e-6, 100e-9]),
    },
    2: {
        "Nfrac": np.array([0.18, 0.74, 1.0e-8]),
        "logSigma": np.array([0.19, 0.45, 0.25]),
        "Dm": np.array([0.122e-6, 0.140e-6, 100e-9]),
    },
    3: {
        "Nfrac": np.array([0.16, 0.91, 1.0e-8]),
        "logSigma": np.array([0.19, 0.43, 0.25]),
        "Dm": np.array([0.084e-6, 0.115e-6, 100e-9]),
    },
    4: {
        "Nfrac": np.array([0.20, 1.06, 1.0e-8]),
        "logSigma": np.array([0.23, 0.47, 0.25]),
        "Dm": np.array([0.061e-6, 0.102e-6, 100e-9]),
    },
    5: {
        "Nfrac": np.array([0.60, 1.37, 1.0e-8]),
        "logSigma": np.array([0.49, 0.76, 0.25]),
        "Dm": np.array([0.038e-6, 0.080e-6, 100e-9]),
    },
    6: {
        "Nfrac": np.array([0.56, 1.118, 1.0e-8]),
        "logSigma": np.array([0.46, 0.68, 0.25]),
        "Dm": np.array([0.029e-6, 0.053e-6, 100e-9]),
    },
}


SCRIPT_DIR = Path(__file__).resolve().parent
NAMELIST = SCRIPT_DIR / "namelist-aida.in"
MODEL_DIR = (SCRIPT_DIR / "../..").resolve()
MODEL_EXE = MODEL_DIR / "main.exe"

USERNAME = getpass.getuser()


def float_tag(value: float) -> str:
    """Filename-safe compact representation of a floating-point value."""
    return f"{value:g}".replace("-", "m").replace(".", "p")


# Keep different updraft experiments separate.  Diagnostics do not need
# separate model runs, so they share the same NetCDF directory.
OUTPUT_DIR = Path("/tmp") / USERNAME / f"aida_psd_sweep_w{float_tag(WINIT)}"


# -----------------------------------------------------------------------------
# Namelist helpers
# -----------------------------------------------------------------------------

def set_namelist_assignment(text: str, name: str, rhs: str) -> str:
    """Replace a one-line Fortran namelist assignment by variable name."""
    pattern = re.compile(
        rf"^(?P<indent>[ \t]*){re.escape(name)}[ \t]*=[^\n]*(?P<newline>\n|$)",
        re.MULTILINE,
    )
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one assignment for {name!r} in the namelist; "
            f"found {len(matches)}."
        )

    def repl(match: re.Match[str]) -> str:
        return f"{match.group('indent')}{name} = {rhs},{match.group('newline')}"

    return pattern.sub(repl, text, count=1)


def replace_namelist_group(text: str, group: str, body: str) -> str:
    """Replace an entire Fortran namelist group, independent of its contents."""
    pattern = re.compile(
        rf"^[ \t]*&{re.escape(group)}\b.*?^[ \t]*/[ \t]*(?:!.*)?$",
        re.MULTILINE | re.DOTALL,
    )
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one &{group} namelist group; found {len(matches)}."
        )
    replacement = f"&{group}\n{body.rstrip()}\n/"
    return pattern.sub(lambda _: replacement, text, count=1)


def format_values(values: np.ndarray) -> str:
    return ", ".join(f"{float(value):.12g}" for value in values)


def mode_numbers_per_kg(
    total_number_m3: float,
    psd_number_fractions: np.ndarray,
) -> np.ndarray:
    """Convert prescribed total number per m3 air to model number per kg air."""
    fractions = psd_number_fractions / np.sum(psd_number_fractions)
    return total_number_m3 * fractions / RHO_AIR_INITIAL


def output_file(psd_type: int, run_index: int) -> Path:
    return OUTPUT_DIR / f"output_psd{psd_type:02d}_{run_index:03d}.nc"


def make_namelist(
    template: str,
    psd_type: int,
    added_number_cm3: float,
    run_index: int,
) -> str:
    """Create one warm, adiabatic-parcel model namelist."""
    background = PSD[BACKGROUND_PSD]
    added = PSD[psd_type]

    n1 = mode_numbers_per_kg(BACKGROUND_NUMBER_CM3 * 1.0e6, background["Nfrac"])
    n2 = mode_numbers_per_kg(added_number_cm3 * 1.0e6, added["Nfrac"])

    text = template

    # This is explicitly an adiabatic parcel experiment, not a chamber run.
    # Keep only the old/common chamber switch so newer chamber options cannot
    # make an older executable reject the namelist.
    text = replace_namelist_group(
        text,
        "chamber_options",
        "    n_levels_c = 0,",
    )

    text = set_namelist_assignment(
        text, "outputfile", f"'{output_file(psd_type, run_index)}'"
    )
    text = set_namelist_assignment(text, "winit", f"{WINIT:.12g}")
    text = set_namelist_assignment(text, "updraft_type", "1")
    text = set_namelist_assignment(text, "adiabatic_prof", ".true.")
    text = set_namelist_assignment(text, "entrain_period", "0")
    text = set_namelist_assignment(text, "ent_rate", "0.0")
    text = set_namelist_assignment(text, "ice_flag", "0")
    text = set_namelist_assignment(text, "sce_flag", "0")
    text = set_namelist_assignment(text, "use_prof_for_tprh", ".false.")
    text = set_namelist_assignment(text, "z_ctop", "-1.0")

    # Background ammonium sulphate: always PSD2.
    text = set_namelist_assignment(
        text, "n_aer1(1:3,1:1)", format_values(n1)
    )
    text = set_namelist_assignment(
        text, "d_aer1(1:3,1:1)", format_values(background["Dm"])
    )
    text = set_namelist_assignment(
        text, "sig_aer1(1:3,1:1)", format_values(background["logSigma"])
    )

    # Added NaCl: all three internal submodes.
    text = set_namelist_assignment(
        text, "n_aer1(1:3,2:2)", format_values(n2)
    )
    text = set_namelist_assignment(
        text, "d_aer1(1:3,2:2)", format_values(added["Dm"])
    )
    text = set_namelist_assignment(
        text, "sig_aer1(1:3,2:2)", format_values(added["logSigma"])
    )

    return text


# -----------------------------------------------------------------------------
# Model execution
# -----------------------------------------------------------------------------

def prescribed_added_numbers() -> np.ndarray:
    return np.linspace(0.0, MAX_ADDED_NUMBER_CM3, N_RUNS)


def run_all(namelist: Path = NAMELIST) -> np.ndarray:
    """Run all six PSD sweeps and return prescribed NaCl concentrations."""
    added_number_cm3 = prescribed_added_numbers()

    if not BMM_RUN:
        return added_number_cm3

    if not namelist.exists():
        raise FileNotFoundError(f"Namelist not found: {namelist}")
    if not MODEL_EXE.exists():
        raise FileNotFoundError(f"Model executable not found: {MODEL_EXE}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    template = namelist.read_text()

    total = len(ADDED_PSDS) * N_RUNS
    count = 0

    for psd_type in ADDED_PSDS:
        print(f"\n=== Added NaCl PSD{psd_type}; background PSD{BACKGROUND_PSD} ===")
        for k, number_cm3 in enumerate(added_number_cm3):
            count += 1
            print(
                f"[{count:3d}/{total}] PSD{psd_type}, run {k:03d}: "
                f"NaCl = {number_cm3:.2f} cm^-3"
            )

            namelist_text = make_namelist(template, psd_type, number_cm3, k)

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".in", delete=False
            ) as tmp:
                tmp.write(namelist_text)
                tmp_name = Path(tmp.name)

            try:
                subprocess.run(
                    [str(MODEL_EXE), str(tmp_name)],
                    cwd=MODEL_DIR,
                    check=True,
                )
            except subprocess.CalledProcessError:
                failed = OUTPUT_DIR / (
                    f"failed_namelist_psd{psd_type:02d}_{k:03d}.in"
                )
                failed.write_text(namelist_text)
                print(f"\nFailed namelist saved to: {failed}")
                raise
            finally:
                tmp_name.unlink(missing_ok=True)

    return added_number_cm3


# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------

def air_density(nc: Dataset) -> np.ndarray:
    return np.asarray(nc["p"][:], dtype=float) / np.asarray(
        nc["t"][:], dtype=float
    ) / R_D


def concentration_cm3(number_per_kg: np.ndarray, rho_air: np.ndarray) -> np.ndarray:
    """Convert a time series of number kg^-1 to number cm^-3."""
    return np.asarray(number_per_kg, dtype=float) * rho_air / 1.0e6


def ndrop_timeseries_cm3(nc: Dataset) -> np.ndarray:
    """Activated droplet number from the model ``ndrop`` field [cm^-3]."""
    ndrop = np.ma.asarray(nc["ndrop"][:])

    # ndrop is normally a 1-D total.  If a future output stores components in
    # trailing dimensions, sum those components at each time rather than
    # averaging them separately.
    if ndrop.ndim > 1:
        ndrop = np.ma.sum(ndrop, axis=tuple(range(1, ndrop.ndim)))

    ndrop = np.ma.filled(ndrop, np.nan).astype(float)
    return concentration_cm3(ndrop, air_density(nc))


def dgt_timeseries_cm3(nc: Dataset, dmin_um: float) -> np.ndarray:
    """Number concentration of liquid drops with D > dmin_um [cm^-3].

    This follows the original lutAIDA.py diagnostic: ``mwat`` is converted to
    an equivalent liquid-water diameter using rho_w=1000 kg m^-3, and ``nwat``
    is summed only where the corresponding diameter exceeds the threshold.
    """
    nwat = np.ma.asarray(nc["nwat"][:])
    mwat = np.ma.asarray(nc["mwat"][:])

    if nwat.shape != mwat.shape:
        raise RuntimeError(
            f"nwat and mwat shapes differ: {nwat.shape} versus {mwat.shape}"
        )
    if nwat.ndim < 2:
        raise RuntimeError(
            f"Expected nwat/mwat to have time plus bin dimensions; got {nwat.shape}"
        )

    # m = rho * pi/6 * D^3
    diameter_m = (mwat / (np.pi / 6.0 * RHO_WATER)) ** (1.0 / 3.0)
    threshold_m = dmin_um * 1.0e-6

    selected = np.ma.where(diameter_m > threshold_m, nwat, 0.0)
    number_per_kg = np.ma.sum(selected, axis=tuple(range(1, selected.ndim)))
    number_per_kg = np.ma.filled(number_per_kg, np.nan).astype(float)

    return concentration_cm3(number_per_kg, air_density(nc))


def drop_number_from_file(
    filename: Path,
    diagnostic: str,
    dmin_um: float,
) -> float:
    """Return time-mean drop number concentration for the selected diagnostic."""
    with Dataset(filename) as nc:
        if diagnostic == "ndrop":
            conc = ndrop_timeseries_cm3(nc)
        elif diagnostic == "dgt":
            conc = dgt_timeseries_cm3(nc, dmin_um)
        else:
            raise ValueError(f"Unknown diagnostic {diagnostic!r}")

        return float(np.nanmean(conc))


def diagnostic_description(diagnostic: str, dmin_um: float) -> str:
    if diagnostic == "ndrop":
        return "activated drops (ndrop)"
    return rf"drops with D > {dmin_um:g} $\mu$m"


def diagnostic_file_tag(diagnostic: str, dmin_um: float) -> str:
    if diagnostic == "ndrop":
        return "ndrop"
    return f"dgt_{float_tag(dmin_um)}um"


def analyse(
    diagnostic: str,
    dmin_um: float,
) -> tuple[np.ndarray, dict[int, np.ndarray]]:
    """Read all outputs and return one drop-number curve per added PSD."""
    added_number_cm3 = prescribed_added_numbers()
    results: dict[int, np.ndarray] = {}

    for psd_type in ADDED_PSDS:
        drop_number = np.zeros(N_RUNS)
        for k in range(N_RUNS):
            filename = output_file(psd_type, k)
            if not filename.exists():
                raise FileNotFoundError(
                    f"Missing model output: {filename}\n"
                    "Run without --plot-only first, or check the failed model run."
                )
            drop_number[k] = drop_number_from_file(
                filename, diagnostic, dmin_um
            )
        results[psd_type] = drop_number

    return added_number_cm3, results


# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

def save_minima(
    added_number_cm3: np.ndarray,
    results: dict[int, np.ndarray],
    diagnostic: str,
    dmin_um: float,
) -> Path:
    """Write the minimum drop number and its NaCl concentration for each PSD."""
    csv_file = OUTPUT_DIR / f"aida_psd_sweep_minima_{diagnostic_file_tag(diagnostic, dmin_um)}.csv"

    with csv_file.open("w", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(
            [
                "added_psd",
                "NaCl_at_minimum_cm-3",
                "minimum_drop_number_cm-3",
                "diagnostic",
                "diameter_threshold_um",
                "background_psd",
                "background_number_cm-3",
                "w_m_s-1",
                "aerosol_number_scale",
            ]
        )
        for psd_type in ADDED_PSDS:
            i_min = int(np.nanargmin(results[psd_type]))
            writer.writerow(
                [
                    psd_type,
                    added_number_cm3[i_min],
                    results[psd_type][i_min],
                    diagnostic,
                    dmin_um if diagnostic == "dgt" else "",
                    BACKGROUND_PSD,
                    BACKGROUND_NUMBER_CM3,
                    WINIT,
                    AEROSOL_NUMBER_SCALE,
                ]
            )

    return csv_file


def make_plot(
    added_number_cm3: np.ndarray,
    results: dict[int, np.ndarray],
    diagnostic: str,
    dmin_um: float,
) -> Path:
    """Plot selected drop-number diagnostic versus prescribed NaCl number."""
    fig, ax = plt.subplots(figsize=(8.0, 5.5))

    for psd_type in ADDED_PSDS:
        y = results[psd_type]
        line, = ax.plot(
            added_number_cm3,
            y,
            marker="o",
            markersize=3,
            linewidth=1.5,
            label=f"NaCl PSD{psd_type}",
        )

        i_min = int(np.nanargmin(y))
        ax.plot(
            added_number_cm3[i_min],
            y[i_min],
            marker="o",
            markersize=7,
            linestyle="none",
            color=line.get_color(),
        )

    ax.set_xlabel(r"Added NaCl aerosol number concentration (cm$^{-3}$)")

    if diagnostic == "ndrop":
        ax.set_ylabel(r"Activated droplet number concentration (cm$^{-3}$)")
    else:
        ax.set_ylabel(
            rf"$N(D>{dmin_um:g}\,\mu\mathrm{{m}})$ (cm$^{{-3}}$)"
        )

    ax.set_title(
        rf"Background ammonium sulphate: PSD{BACKGROUND_PSD}, "
        rf"{BACKGROUND_NUMBER_CM3:.1f} cm$^{{-3}}$; "
        rf"$w={WINIT:g}$ m s$^{{-1}}$"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()

    figure_file = OUTPUT_DIR / f"aida_psd_sweep_{diagnostic_file_tag(diagnostic, dmin_um)}.png"
    fig.savefig(figure_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return figure_file


def print_minima(
    added_number_cm3: np.ndarray,
    results: dict[int, np.ndarray],
    diagnostic: str,
    dmin_um: float,
) -> None:
    print(f"\nMinimum {diagnostic_description(diagnostic, dmin_um)} for each added NaCl PSD")
    print("PSD   NaCl at minimum (cm^-3)   minimum N (cm^-3)")
    print("---   -----------------------   -----------------")
    for psd_type in ADDED_PSDS:
        i_min = int(np.nanargmin(results[psd_type]))
        print(
            f"{psd_type:>3d}   {added_number_cm3[i_min]:>23.2f}   "
            f"{results[psd_type][i_min]:>17.2f}"
        )


def print_experiment_summary(diagnostic: str, dmin_um: float) -> None:
    print("Adiabatic parcel PSD sweep")
    print(f"  w = {WINIT:g} m s^-1")
    if SCALE_AEROSOL_WITH_UPDRAFT:
        print(
            "  aerosol number scaling = "
            f"({WINIT:g}/{REFERENCE_UPDRAFT_MS:g})^"
            f"{AEROSOL_SCALING_EXPONENT:g} = {AEROSOL_NUMBER_SCALE:.6f}"
        )
    else:
        print("  aerosol number scaling = OFF")
    print(f"  background PSD{BACKGROUND_PSD} = {BACKGROUND_NUMBER_CM3:.2f} cm^-3")
    print(f"  added NaCl sweep = 0--{MAX_ADDED_NUMBER_CM3:.2f} cm^-3")
    print(f"  diagnostic = {diagnostic_description(diagnostic, dmin_um)}")
    print(f"  output directory = {OUTPUT_DIR}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run and plot adiabatic-parcel NaCl sweeps for added PSD1--PSD6."
    )
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Do not run the model; only analyse existing NetCDF files and replot.",
    )
    parser.add_argument(
        "--namelist",
        type=Path,
        default=NAMELIST,
        help="Template namelist to use (default: namelist-aida.in beside this script).",
    )
    parser.add_argument(
        "--diagnostic",
        choices=("ndrop", "dgt"),
        default=DROP_DIAGNOSTIC,
        help=(
            "Drop-number diagnostic: 'ndrop' uses the model activated-drop field; "
            "'dgt' counts drops larger than --dmin-um."
        ),
    )
    parser.add_argument(
        "--dmin-um",
        type=float,
        default=DROP_DIAMETER_THRESHOLD_UM,
        help="Diameter threshold in microns for --diagnostic dgt (default: 2).",
    )
    args = parser.parse_args()

    if args.dmin_um <= 0.0:
        parser.error("--dmin-um must be > 0")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print_experiment_summary(args.diagnostic, args.dmin_um)

    if not args.plot_only:
        run_all(args.namelist)

    added_number_cm3, results = analyse(args.diagnostic, args.dmin_um)
    csv_file = save_minima(
        added_number_cm3, results, args.diagnostic, args.dmin_um
    )
    figure_file = make_plot(
        added_number_cm3, results, args.diagnostic, args.dmin_um
    )
    print_minima(added_number_cm3, results, args.diagnostic, args.dmin_um)

    print(f"\nFigure: {figure_file}")
    print(f"Minima table: {csv_file}")


if __name__ == "__main__":
    main()
