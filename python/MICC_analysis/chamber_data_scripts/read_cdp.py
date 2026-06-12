import pandas as pd
import numpy as np

path1='/Volumes/REFLECT/data/20260512/cdp/20260512115456/'
path1='/Volumes/REFLECT/data/20260526/cdp/20260526152153_exp2/'


import re
import pandas as pd


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
        raise ValueError("Could not find data header line starting with 'End Seconds'.")

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

    df["elapsed_seconds"] = (
        df["datetime"] - df["datetime"].iloc[0]
    ).dt.total_seconds()

    return {
        "data": df,
        "sizes_um": sizes_um,
        "thresholds": thresholds,
        "ipt_thresholds": ipt_thresholds,
        "sample_area_mm2": sample_area_mm2,
        "sample_time_s": sample_time_s,
        "metadata": metadata,
    }


filename = path1 + "00CDP PBP20260526152153.csv"

cdp = read_cdp_pbp(filename)


df = cdp["data"]
sizes_um=np.array(cdp['sizes_um'])
bin_widths=np.zeros_like(sizes_um)
ind,=np.where(sizes_um<=14)
bin_widths[ind]=2
ind,=np.where(sizes_um>14)
bin_widths[ind]=3


# conc in bins
cdp_conc2=np.zeros((len(df['End Seconds']),30))
for i in range(30):
	cdp_conc2[:,i]=df['CDP Bin ' + str(i+1)]/ \
		(cdp["sample_area_mm2"]*1e-6*10.0*cdp["sample_time_s"])


lwc=np.pi/6*1000*np.sum(cdp_conc2*(sizes_um*1e-6)**3,axis=1)*1e3

print("Sizes / um:")
print(sizes_um)

print("\nThresholds:")
print(cdp['thresholds'])

print("\nIPT thresholds:")
print(cdp['ipt_thresholds'])

print("\nData:")
print(df.head())