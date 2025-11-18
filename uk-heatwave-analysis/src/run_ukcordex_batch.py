from pathlib import Path
import numpy as np
import pandas as pd

from detect_heatwaves import (
    load_region_mask_indices,
    process_one_pair
)


def main():
    """
    Example driver script.

    It:
        1. Scans the UKCORDEX directory tree for valid model folders
           (similar to the FolderList filtering in the original R code).
        2. Loads the OSGB region mask and defines temperature thresholds.
        3. For each (GCM, RCM, RunID) combination, finds matching tasmax/tasmin
           daily files, runs the heatwave analysis for all regions, and
           aggregates results across periods.
        4. Writes a CSV file 'HWInfo_{GCM}_{RCM}_{RunID}.csv' with annual
           regional heatwave statistics.
    """

    root = Path("/data/met/ukcordex")

    folder_list = []
    for p in root.rglob("*"):
        if not p.is_dir():
            continue
        s = str(p)
        if "ceda" not in s:
            continue
        if "checkpoints" in s:
            continue
        if "ECMWF-ERAINT" in s:
            continue
        folder_list.append(p)

    folder_info = []
    for folder in folder_list:
        parts = folder.parts

        try:
            gcm = parts[4]
            rcm = parts[5]
            runid = parts[6]
        except IndexError:
            continue
        folder_info.append({"folder": folder, "GCM": gcm, "RCM": rcm, "RunID": runid})

    folder_info = pd.DataFrame(folder_info)

    # ---- Load region mask and thresholds ----
    mask_path = "regionmask-region_osgb.nc"
    regions_list = load_region_mask_indices(mask_path)

    thresholds = pd.DataFrame({
        "region": [info["Region_Name"] for _, info in regions_list.items()],
        "max_threshold": [30, 30, 28, 32, 28, 28, 30, 31, 30, 30, 28, 29, 30, 30, 28, 30],
        "min_threshold": [15, 15, 15, 18, 15, 15, 15, 16, 15, 15, 15, 15, 15, 15, 15, 15],
    })  # Threshold is published by PHE(Public Health England)
    # IMPORTANT: ensure the order in these threshold lists matches the order
    # of Region_Name in regions_list. You may want to reorder or construct
    # this table explicitly based on the region names.

    for _, row in folder_info.iterrows():
        folder = row["folder"]
        gcm = row["GCM"]
        rcm = row["RCM"]
        runid = row["RunID"]

        out_name = f"HWInfo_{gcm}_{rcm}_{runid}.csv"
        if Path(out_name).exists():
            print(f"{out_name} already exists, skip.")
            continue

        print(f"Processing data in folder {folder} ...")

        tasmax_files = sorted([f for f in folder.glob("*tasmax*day*12km*")])
        tasmin_files = sorted([f for f in folder.glob("*tasmin*day*12km*")])

        if len(tasmax_files) != len(tasmin_files):
            print("tasmax/tasmin file count mismatch, skip this folder.")
            continue


        all_periods = []
        for tmax_path, tmin_path in zip(tasmax_files, tasmin_files):
            print(f"  Pair: {tmax_path.name} / {tmin_path.name}")
            df_period = process_one_pair(
                tasmax_path=str(tmax_path),
                tasmin_path=str(tmin_path),
                regions_list=regions_list,
                thresholds_df=thresholds,
                min_consecutive_days=3,
            )
            all_periods.append(df_period)

        if not all_periods:
            print("No data generated, skip.")
            continue

        heatwave_summary = pd.concat(all_periods, ignore_index=True)

        # Aggregate across multiple periods for the same year
        def adjust(x):
            x = np.asarray(x)
            if np.any(x == 0):
                return x.sum()
            else:
                return x.mean()

        summary_final = (
            heatwave_summary
            .groupby(["Year", "Region_Name", "Region_ID"])
            .agg({
                "Heatwave_Count": adjust,
                "Mean_Duration": adjust,
                "Mean_Intensity": adjust,
            })
            .reset_index()
        )

        print(f"Writing to {out_name}")
        summary_final.to_csv(out_name, index=False)


if __name__ == "__main__":
    main()
