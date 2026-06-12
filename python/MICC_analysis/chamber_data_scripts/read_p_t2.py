import re
from pathlib import Path

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import svp

r_gas=8.314
mair=29e-3
mh2o=18e-3
eps1=mh2o/mair
rair=r_gas/mair

def add_datetime_and_seconds_in_day(df):
    df.columns = df.columns.str.strip()

    df["datetime"] = pd.to_datetime(
        df["Date"].astype(str).str.strip()
        + " "
        + df["Time"].astype(str).str.strip(),
        format="%d/%m/%y %H:%M:%S"
    )

    df["seconds_in_day"] = (
        df["datetime"].dt.hour * 3600
        + df["datetime"].dt.minute * 60
        + df["datetime"].dt.second
        + df["datetime"].dt.microsecond / 1e6
    )

    return df


def read_thermocouple_csv(filename):
    df = pd.read_csv(filename, skiprows=1)

    # Remove empty columns caused by trailing commas
    df = df.loc[:, ~df.columns.str.contains("^Unnamed")]

    df = add_datetime_and_seconds_in_day(df)

    return df


def read_keller_s30(filename):
    df = pd.read_csv(filename, sep="\t", skiprows=1)

    df = add_datetime_and_seconds_in_day(df)

    df["Pressure (mBar)"] = pd.to_numeric(
        df["Pressure (mBar)"],
        errors="coerce"
    )

    return df


def find_temperature_pressure_files(base_path, date_string=None):
    """
    Finds files of the form:
        YYYYMMDD-temp.csv
        YYYYMMDD-pressure.S30
    """

    base_path = Path(base_path)

    temp_pattern = re.compile(r"^(\d{8})-temp\.csv$")
    pressure_pattern = re.compile(r"^(\d{8})-pressure\.S30$")

    temp_files = []
    pressure_files = []

    for f in base_path.iterdir():
        if not f.is_file():
            continue

        mt = temp_pattern.match(f.name)
        mp = pressure_pattern.match(f.name)

        if mt:
            if date_string is None or mt.group(1) == date_string:
                temp_files.append(f)

        if mp:
            if date_string is None or mp.group(1) == date_string:
                pressure_files.append(f)

    temp_files = sorted(temp_files)
    pressure_files = sorted(pressure_files)

    if date_string is None:
        temp_dates = {
            temp_pattern.match(f.name).group(1)
            for f in temp_files
        }

        pressure_dates = {
            pressure_pattern.match(f.name).group(1)
            for f in pressure_files
        }

        common_dates = sorted(temp_dates & pressure_dates)

        if len(common_dates) == 0:
            raise FileNotFoundError(
                f"No matching YYYYMMDD-temp.csv and "
                f"YYYYMMDD-pressure.S30 pair found in {base_path}"
            )

        if len(common_dates) > 1:
            raise ValueError(
                "More than one matching date found: "
                + ", ".join(common_dates)
                + ". Pass date_string='YYYYMMDD' to choose one."
            )

        date_string = common_dates[0]

        temp_files = [
            f for f in temp_files
            if f.name.startswith(date_string)
        ]

        pressure_files = [
            f for f in pressure_files
            if f.name.startswith(date_string)
        ]

    if len(temp_files) != 1:
        raise FileNotFoundError(
            f"Expected one temperature file for {date_string}, "
            f"found {len(temp_files)}"
        )

    if len(pressure_files) != 1:
        raise FileNotFoundError(
            f"Expected one pressure file for {date_string}, "
            f"found {len(pressure_files)}"
        )

    return temp_files[0], pressure_files[0], date_string


def interpolate_pressure_to_thermo_time(thermo, pressure):
    """
    Interpolates pressure onto the thermocouple timebase.

    Converts pressure from mBar to Pa.
    """

    valid = np.isfinite(pressure["Pressure (mBar)"])

    intp = interp1d(
        np.array(pressure.loc[valid, "seconds_in_day"]),
        np.array(pressure.loc[valid, "Pressure (mBar)"]) * 100.0,
        bounds_error=False,
        fill_value="extrapolate"
    )

    press_temp = intp(thermo["seconds_in_day"])

    return press_temp, intp


def array_to_fortran_line(name, arr):
    values = ",".join(f"{x:.8g}" for x in arr)
    return f"{name}(1:{len(arr)})={values},\n"


def make_chamber_output(thermo, pressure, t0, rhinit, tlen):
    """
    Makes chamber data and Fortran-style string output.
    """

    press_temp, intp = interpolate_pressure_to_thermo_time(
        thermo,
        pressure
    )


    # Average T3, T4, T5 and convert from degC to K
    temperature = (
        thermo["T3 (C)"]
        + thermo["T4 (C)"]
        + thermo["T5 (C)"]
    ) / 3 + 273.15

    # Interpolate temperature relative to t0
    intt = interp1d(
        np.array(thermo["seconds_in_day"] - t0),
        np.array(temperature),
        bounds_error=False,
        fill_value="extrapolate"
    )

    tinit = intt(0)
    pinit = intp(t0)

    es = svp.svp([tinit], "buck2", "liq")[0]

    qtot0 = rhinit * es / (pinit - es) * eps1
    qtot = qtot0 * np.ones(len(temperature))

    inds, = np.where(
        (np.array(thermo["seconds_in_day"] - t0) >= 0)
        & (np.array(thermo["seconds_in_day"] - t0) <= tlen)
    )

    time_chamber = np.array(
        thermo["seconds_in_day"].iloc[inds] - t0
    )

    press_chamber = np.array(press_temp[inds])
    temp_chamber = np.array(temperature.iloc[inds])
    qtot_chamber = np.array(qtot[inds])

    str1 = ""
    str1 += array_to_fortran_line("time_chamber", time_chamber)
    str1 += array_to_fortran_line("press_chamber", press_chamber)
    str1 += array_to_fortran_line("temp_chamber", temp_chamber)
    str1 += array_to_fortran_line("qtot_chamber", qtot_chamber)

    chamber = {
        "t0": t0,
        "tlen": tlen,
        "rhinit": rhinit,
        "tinit": tinit,
        "pinit": pinit,
        "es": es,
        "qtot0": qtot0,
        "time_chamber": time_chamber,
        "press_chamber": press_chamber,
        "temp_chamber": temp_chamber,
        "qtot_chamber": qtot_chamber,
    }

    return str1, chamber, press_temp, intp, temperature, qtot


def read_temperature_pressure(
    base_path,
    date_string=None,
    t0=None,
    rhinit=None,
    tlen=None
):
    """
    Reads temperature and pressure files from a directory.

    Files should be named:
        YYYYMMDD-temp.csv
        YYYYMMDD-pressure.S30

    If t0, rhinit, and tlen are provided, also creates str1.

    Returns
    -------
    tp : dict
        Dictionary containing all data and metadata.
    str1 : str
        Fortran-style chamber string.
    """

    temp_file, pressure_file, date_string = find_temperature_pressure_files(
        base_path,
        date_string=date_string
    )

    thermo = read_thermocouple_csv(temp_file)
    pressure = read_keller_s30(pressure_file)

    tp = {
        "thermo": thermo,
        "pressure": pressure,
        "temp_file": temp_file,
        "pressure_file": pressure_file,
        "date_string": date_string,
        "str1": "",
    }

    str1 = ""

    if t0 is not None and rhinit is not None and tlen is not None:
        (
            str1,
            chamber,
            press_temp,
            pressure_interpolator,
            temperature,
            qtot,
        ) = make_chamber_output(
            thermo=thermo,
            pressure=pressure,
            t0=t0,
            rhinit=rhinit,
            tlen=tlen
        )

        tp["str1"] = str1
        tp["chamber"] = chamber
        tp["press_temp"] = press_temp
        tp["pressure_interpolator"] = pressure_interpolator
        tp["temperature"] = temperature
        tp["qtot"] = qtot

    return tp, str1