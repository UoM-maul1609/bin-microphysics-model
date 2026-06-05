#!/usr/bin/env python3
"""
Read a UHSAS tab-delimited "xls" export where:
- Row 1 contains column labels, with bin LOWER edges as numeric strings (e.g. 60.00, 61.73, 63.51, ...)
- Row 2 contains units for metadata columns, with bin UPPER edges as numeric strings (e.g. 61.73, 63.51, ...)

Outputs:
- df_meta: metadata dataframe
- counts: 2D numpy array (n_times, n_bins)
- bin_edges: 1D numpy array (n_bins+1)
- bin_lower, bin_upper: 1D arrays (n_bins)
- time: pandas datetime Series
"""

from __future__ import annotations

import csv
from pathlib import Path
import numpy as np
import pandas as pd


def _first_numeric_col(tokens: list[str]) -> int:
    """Return index of first token that can be parsed as float; raise if none."""
    for i, t in enumerate(tokens):
        try:
            float(t)
            return i
        except Exception:
            pass
    raise ValueError("Could not find first numeric (bin-edge) column in header row.")


def read_uhsas_file(file_path: str | Path):
    file_path = Path(file_path)

    # Read first two header rows (tab-delimited)
    with file_path.open("r", newline="", errors="replace") as f:
        reader = csv.reader(f, delimiter="\t")
        header1 = next(reader)
        header2 = next(reader)

    bin_start = _first_numeric_col(header1)

    meta_cols = header1[:bin_start]
    bin_lower = np.array([float(x) for x in header1[bin_start:]], dtype=float)
    bin_upper = np.array([float(x) for x in header2[bin_start:]], dtype=float)

    # Build bin edges (Nbins+1): [first lower, then all uppers]
    bin_edges = np.concatenate(([bin_lower[0]], bin_upper))

    # Read the data body
    # Use names from header1 so bin columns are named by their LOWER edge strings
    df = pd.read_csv(
        file_path,
        sep="\t",
        skiprows=2,
        names=header1,
        engine="python",
    )

    # Parse datetime from Date + Time columns if present
    if "Date" in df.columns and "Time" in df.columns:
        # Some exports have double spaces before AM/PM; strip whitespace.
        dt_str = (df["Date"].astype(str).str.strip() + " " + df["Time"].astype(str).str.strip())
        # Example: 2/24/2026 04:49:10.29 AM  (fractional seconds with 2 decimals)
        # Let pandas infer robustly.
        time = pd.to_datetime(dt_str, errors="coerce")
    else:
        time = pd.Series(pd.NaT, index=df.index)

    # Split metadata vs counts
    df_meta = df[meta_cols].copy() if meta_cols else pd.DataFrame(index=df.index)

    # Bin columns are everything from bin_start onward (named by header1)
    bin_col_names = header1[bin_start:]
    counts = df[bin_col_names].apply(pd.to_numeric, errors="coerce").to_numpy()

    return {
        "df_meta": df_meta,
        "time": time,
        "counts": counts,
        "bin_edges": bin_edges,
        "bin_lower": bin_lower,
        "bin_upper": bin_upper,
        "bin_col_names": bin_col_names,
    }


if __name__ == "__main__":
    # Change this to your file name/path
    path = "20260224044900.xls"

    out = read_uhsas_file(path)

    print("Rows:", out["counts"].shape[0])
    print("Bins:", out["counts"].shape[1])
    print("First 5 bin edges:", out["bin_edges"][:5])
    print("Last 5 bin edges:", out["bin_edges"][-5:])
    print("\nMetadata columns:", list(out["df_meta"].columns))
    print("First timestamps:\n", out["time"].head())
    