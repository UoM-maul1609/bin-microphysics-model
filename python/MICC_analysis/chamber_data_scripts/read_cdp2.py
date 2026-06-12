import re
from pathlib import Path

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d


def parse_array_from_metadata(line):
    values = line.split("=", 1)[1].strip()
    values = re.sub(r"^<\d+>", "", values)
    return [float(x.strip()) for x in values.split(",") if x.strip() != ""]


def parse_numeric_metadata_value(value):
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", value)
    if match is None:
        return None
    return float(match.group(0))


def read_cdp_pbp(filename):
    metadata = {}

    sizes_um = None
    thresholds = None
    ipt_thresholds = None
    sample_area_mm2 = None
    sample_time_s = None
    header_line = None

    filename = Path(filename)

    with open(filename, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        line_clean = line.strip()

        if line_clean.startswith("End Seconds"):
            header_line = i
            break

        if "=" in line_clean:
            key, value = line_clean.split("=", 1)
            key = key.strip()
            value = value.strip()

            metadata[key] = value

            if key == "Sizes":
                sizes_um = parse_array_from_metadata(line_clean)

            elif key == "Thresholds":
                thresholds = parse_array_from_metadata(line_clean)

            elif key == "IPT Thresholds":
                ipt_thresholds = parse_array_from_metadata(line_clean)

            elif key.startswith("Sample Area"):
                sample_area_mm2 = parse_numeric_metadata_value(value)

            elif key.startswith("Sample Time"):
                sample_time_s = parse_numeric_metadata_value(value)

    if header_line is None:
        raise ValueError(f"Could not find data header line in {filename}")

    df = pd.read_csv(filename, skiprows=header_line)
    df.columns = df.columns.str.strip()

    df["Year"] = pd.to_numeric(df["Year"], errors="coerce")
    df["Day of Year"] = pd.to_numeric(df["Day of Year"], errors="coerce")
    df["End Seconds"] = pd.to_numeric(df["End Seconds"], errors="coerce")

    df["seconds_in_day"] = df["End Seconds"]

    df["datetime"] = (
        pd.to_datetime(df["Year"].astype("Int64").astype(str), format="%Y")
        + pd.to_timedelta(df["Day of Year"] - 1, unit="D")
        + pd.to_timedelta(df["seconds_in_day"], unit="s")
    )

    df["source_file"] = str(filename)

    return {
        "data": df,
        "sizes_um": sizes_um,
        "thresholds": thresholds,
        "ipt_thresholds": ipt_thresholds,
        "sample_area_mm2": sample_area_mm2,
        "sample_time_s": sample_time_s,
        "metadata": metadata,
    }


def find_cdp_pbp_files(base_path):
    """
    Recursively find files matching exactly:
        00CDP PBPYYYYMMDDHHMMSS.csv

    Excludes:
        00CDP PBPYYYYMMDDHHMMSS_PBP.csv
    """

    base_path = Path(base_path)

    pattern = re.compile(r"^00CDP PBP\d{14}\.csv$")

    files = [
        f for f in base_path.rglob("00CDP PBP*.csv")
        if pattern.match(f.name)
    ]

    # Sort by timestamp in filename
    files = sorted(files, key=lambda f: f.name)

    return files


def make_nan_gap_row(df):
    """
    Make a one-row dataframe with the same columns as df,
    filled with NaN/NaT, to force plotting gaps between files.
    """

    gap = pd.DataFrame([{col: np.nan for col in df.columns}])

    # Keep datetime as NaT rather than float NaN
    if "datetime" in gap.columns:
        gap["datetime"] = pd.NaT

    # Useful label so you can identify gap rows later
    if "source_file" in gap.columns:
        gap["source_file"] = "NaN gap"

    return gap


def read_all_cdp_pbp(base_path, add_nan_gaps=True):
    files = find_cdp_pbp_files(base_path)

    if len(files) == 0:
        raise FileNotFoundError(
            f"No CDP PBP files found below {base_path}"
        )

    all_dfs = []

    sizes_um = None
    thresholds = None
    ipt_thresholds = None
    sample_area_mm2 = None
    sample_time_s = None

    metadata_by_file = {}

    for file_number, file in enumerate(files):
        cdp = read_cdp_pbp(file)

        df = cdp["data"].copy()
        all_dfs.append(df)

        # Add a NaN row between files, but not after the final file
        if add_nan_gaps and file_number < len(files) - 1:
            all_dfs.append(make_nan_gap_row(df))

        metadata_by_file[str(file)] = cdp["metadata"]

        # Use first file as reference
        if sizes_um is None:
            sizes_um = np.array(cdp["sizes_um"])
            thresholds = np.array(cdp["thresholds"])
            ipt_thresholds = np.array(cdp["ipt_thresholds"])
            sample_area_mm2 = cdp["sample_area_mm2"]
            sample_time_s = cdp["sample_time_s"]

        else:
            # Check metadata consistency
            if not np.allclose(sizes_um, np.array(cdp["sizes_um"])):
                print(f"WARNING: Sizes differ in {file}")

            if not np.allclose(thresholds, np.array(cdp["thresholds"])):
                print(f"WARNING: Thresholds differ in {file}")

            if not np.allclose(ipt_thresholds, np.array(cdp["ipt_thresholds"])):
                print(f"WARNING: IPT thresholds differ in {file}")

            if sample_area_mm2 != cdp["sample_area_mm2"]:
                print(
                    f"WARNING: Sample area differs in {file}: "
                    f"{cdp['sample_area_mm2']} vs {sample_area_mm2}"
                )

            if sample_time_s != cdp["sample_time_s"]:
                print(
                    f"WARNING: Sample time differs in {file}: "
                    f"{cdp['sample_time_s']} vs {sample_time_s}"
                )

    # Concatenate as one long time series
    df_all = pd.concat(all_dfs, ignore_index=True)

    # Do NOT sort after inserting NaN gaps, otherwise the NaT rows move
    # Instead, files are already sorted by filename/time.

    # Recalculate elapsed seconds across the whole combined series.
    # Gap rows stay NaN.
    valid_time = df_all["datetime"].notna()

    df_all["elapsed_seconds"] = np.nan
    df_all.loc[valid_time, "elapsed_seconds"] = (
        df_all.loc[valid_time, "datetime"]
        - df_all.loc[valid_time, "datetime"].iloc[0]
    ).dt.total_seconds()

    return {
        "data": df_all,
        "files": files,
        "sizes_um": sizes_um,
        "thresholds": thresholds,
        "ipt_thresholds": ipt_thresholds,
        "sample_area_mm2": sample_area_mm2,
        "sample_time_s": sample_time_s,
        "metadata_by_file": metadata_by_file,
    }

def add_air_density_to_cdp_timebase(cdp, tp, R_d=287.05):
    """
    Interpolates chamber temperature and pressure onto the CDP time-base
    and calculates dry-air density.

    Inputs
    ------
    cdp : dict
        Output from read_all_cdp_pbp().
        Must contain cdp["data"] with seconds_in_day.
    tp : dict
        Output from read_temperature_pressure().
        Must contain:
            tp["thermo"]
            tp["temperature"]
            tp["press_temp"]

    Returns
    -------
    cdp : dict
        Same CDP dictionary, with extra columns added to cdp["data"]:
            pressure_Pa
            temperature_K
            rho_air_kg_m3
    """

    df_cdp = cdp["data"]
    thermo = tp["thermo"]

    # Valid thermo rows
    valid_t = (
        np.isfinite(thermo["seconds_in_day"])
        & np.isfinite(tp["temperature"])
        & np.isfinite(tp["press_temp"])
    )

    # Interpolators from thermo time-base to CDP time-base
    int_temperature = interp1d(
        np.array(thermo.loc[valid_t, "seconds_in_day"]),
        np.array(tp["temperature"][valid_t]),
        bounds_error=False,
        fill_value="extrapolate"
    )

    int_pressure = interp1d(
        np.array(thermo.loc[valid_t, "seconds_in_day"]),
        np.array(tp["press_temp"][valid_t]),
        bounds_error=False,
        fill_value="extrapolate"
    )

    # Only interpolate for real CDP rows, not artificial NaN gap rows
    valid_cdp = np.isfinite(df_cdp["seconds_in_day"])

    df_cdp["temperature_K"] = np.nan
    df_cdp["pressure_Pa"] = np.nan
    df_cdp["rho_air_kg_m3"] = np.nan

    df_cdp.loc[valid_cdp, "temperature_K"] = int_temperature(
        df_cdp.loc[valid_cdp, "seconds_in_day"]
    )

    df_cdp.loc[valid_cdp, "pressure_Pa"] = int_pressure(
        df_cdp.loc[valid_cdp, "seconds_in_day"]
    )

    #R_d = 287.05  # J kg-1 K-1

    df_cdp.loc[valid_cdp, "rho_air_kg_m3"] = (
        df_cdp.loc[valid_cdp, "pressure_Pa"]
        / (R_d * df_cdp.loc[valid_cdp, "temperature_K"])
    )

    cdp["data"] = df_cdp

    return cdp

def add_cdp_concentration_and_lwc(cdp, airspeed_m_s=10.0, rho_water=1000.0):
    """
    Calculates CDP bin concentrations and liquid water content.

    Adds to cdp:
        cdp["bin_widths_um"]
        cdp["sample_volume_m3"]
        cdp["cdp_conc_m3"]
        cdp["cdp_conc_cm3"]

    Adds to cdp["data"]:
        CDP Conc Bin 1 (#/m3), ...
        CDP Conc Bin 1 (#/cm3), ...
        LWC_g_m3

    Inputs
    ------
    cdp : dict
        Output from read_all_cdp_pbp().
    airspeed_m_s : float
        Air speed through the CDP sample area, m/s.
    rho_water : float
        Liquid water density, kg/m3.

    Returns
    -------
    cdp : dict
        Updated CDP dictionary.
    """

    df = cdp["data"]
    sizes_um = np.array(cdp["sizes_um"], dtype=float)

    n_bins = len(sizes_um)

    # Bin widths, matching your previous approach
    bin_widths_um = np.zeros_like(sizes_um)

    ind, = np.where(sizes_um <= 14)
    bin_widths_um[ind] = 2

    ind, = np.where(sizes_um > 14)
    bin_widths_um[ind] = 3

    # Sample volume per CDP sample
    sample_area_m2 = cdp["sample_area_mm2"] * 1e-6
    sample_time_s = cdp["sample_time_s"]

    sample_volume_m3 = sample_area_m2 * airspeed_m_s * sample_time_s
    sample_volume_cm3 = sample_volume_m3 * 1e6

    # Arrays for concentrations
    cdp_conc_m3 = np.full((len(df), n_bins), np.nan)
    cdp_conc_cm3 = np.full((len(df), n_bins), np.nan)

    # Valid rows: avoids artificial NaN gap rows
    valid_rows = np.isfinite(df["seconds_in_day"])

    for i in range(n_bins):
        bin_name = f"CDP Bin {i + 1}"

        counts = pd.to_numeric(df[bin_name], errors="coerce")

        cdp_conc_m3[valid_rows, i] = (
            counts.loc[valid_rows] / sample_volume_m3
        )

        cdp_conc_cm3[valid_rows, i] = (
            counts.loc[valid_rows] / sample_volume_cm3
        )

        df[f"CDP Conc Bin {i + 1} (#/m3)"] = cdp_conc_m3[:, i]
        df[f"CDP Conc Bin {i + 1} (#/cm3)"] = cdp_conc_cm3[:, i]

    # LWC in g/m3
    # cdp_conc_m3 is #/m3
    # diameter is m
    # rho_water * pi/6 * D^3 gives kg per droplet
    lwc_g_m3 = (
        np.pi / 6
        * rho_water
        * np.nansum(cdp_conc_m3 * (sizes_um * 1e-6) ** 3, axis=1)
        * 1e3
    )

    # Keep artificial gap rows as NaN, not zero
    lwc_g_m3[~valid_rows] = np.nan

    df["LWC_g_m3"] = lwc_g_m3

    cdp["data"] = df
    cdp["bin_widths_um"] = bin_widths_um
    cdp["sample_volume_m3"] = sample_volume_m3
    cdp["sample_volume_cm3"] = sample_volume_cm3
    cdp["cdp_conc_m3"] = cdp_conc_m3
    cdp["cdp_conc_cm3"] = cdp_conc_cm3
    cdp["airspeed_m_s"] = airspeed_m_s

    return cdp