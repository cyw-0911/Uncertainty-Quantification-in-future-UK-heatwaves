import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path


def detect_heatwaves_1d(years, max_temp, min_temp,
                        max_threshold, min_threshold,
                        min_consecutive_days = 3):
    """
    Detect heatwaves for a single grid cell (1D time series).

    This function follows the logic of Heatwave Definition:
        1) A heatwave occurs when BOTH max_temp and min_temp exceed
           their respective thresholds.
        2) A valid heatwave must last for at least `min_consecutive_days`.
        3) For each detected heatwave segment, intensity is defined using:
              - the day with the highest max_temp (peak), and
              - the mean max_temp during the segment.
           Intensity = 0.5 * (mean_max_temp + peak_max_temp)

    Parameters
    ----------
    years : Year corresponding to each day in the time series.
    max_temp : Daily maximum temperatures.
    min_temp : Daily minimum temperatures.
    max_threshold : Threshold for maximum temperature.
    min_threshold : Threshold for minimum temperature.
    min_consecutive_days : int, optional
        Minimum duration required for a heatwave event.

    Returns:
    pandas.DataFrame
        Contains, for each year:
            - Year
            - Heatwave_Count
            - Mean_Duration
            - Mean_Intensity
        If no heatwaves are found, returns an empty dataframe.
    """

    max_temp = np.asarray(max_temp)
    min_temp = np.asarray(min_temp)
    years = np.asarray(years)

    # Heatwave condition: both max and min temp thresholds are exceeded
    cond = (max_temp >= max_threshold) & (min_temp >= min_threshold)

    if cond.sum() == 0:
        # No heatwaves at all
        return pd.DataFrame(columns=["Year", "Heatwave_Count", "Mean_Duration", "Mean_Intensity"])

    values = cond.astype(int)
    n = len(values)

    # Indices where the state changes
    change_idx = np.where(np.diff(values) != 0)[0]
    run_ends = np.append(change_idx, n - 1)
    run_starts = np.append(0, run_ends[:-1] + 1)
    run_lengths = run_ends - run_starts + 1
    run_values = values[run_ends]

    events = []

    for start, end, length, val in zip(run_starts, run_ends, run_lengths, run_values):
        # Only keep segments that meet the heatwave criteria
        if val == 1 and length >= min_consecutive_days:
            idx = slice(start, end + 1)
            segment_years = years[idx]
            segment_max = max_temp[idx]
            segment_min = min_temp[idx]

            duration = length

            max_idx_rel = np.argmax(segment_max)
            max_year = segment_years[max_idx_rel]

            mean_max_temp = segment_max.mean()
            peak_max_temp = segment_max[max_idx_rel]

            # Intensity metric Definition 
            intensity = 0.5 * (mean_max_temp + peak_max_temp)

            events.append(
                {
                    "Year": max_year,
                    "Duration": duration,
                    "Intensity": intensity,
                }
            )

    if not events:
        # No segments satisfy the heatwave criteria
        return pd.DataFrame(columns=["Year", "Heatwave_Count", "Mean_Duration", "Mean_Intensity"])

    events_df = pd.DataFrame(events)


    grouped = events_df.groupby("Year").agg(
        Heatwave_Count=("Duration", "count"),
        Total_Duration=("Duration", "sum"),
        Total_Intensity=("Intensity", "sum"),
    ).reset_index()

    grouped["Mean_Duration"] = grouped["Total_Duration"] / grouped["Heatwave_Count"]
    grouped["Mean_Intensity"] = grouped["Total_Intensity"] / grouped["Heatwave_Count"]

    return grouped[["Year", "Heatwave_Count", "Mean_Duration", "Mean_Intensity"]]



def analyze_region_from_arrays(years, tasmax_region, tasmin_region,
                               region_name, region_id,
                               max_threshold, min_threshold,
                               min_consecutive_days=3):
    """
    Compute heatwave metrics for a single region that contains multiple grid cells.

    Parameters
    years : Year
    tasmax_region : Daily maximum temperature for all grid cells in the region.
    tasmin_region : Daily minimum temperature for all grid cells in the region.
    region_name : str, Name of the region.
    region_id : int, Identifier of the region.
    max_threshold : Threshold for maximum temperature.
    min_threshold : Threshold for minimum temperature.
    min_consecutive_days : int, Minimum length (in days) of a heatwave.

    Returns:
    pandas.DataFrame
        One row per year with the following columns:
            - Year
            - Heatwave_Count  (can be averaged over grid cells)
            - Mean_Duration
            - Mean_Intensity
            - Region_Name
            - Region_ID
        If no heatwaves are detected in the region, an empty DataFrame with the same columns is returned.
    """

    years = np.asarray(years)
    n_grid, n_time = tasmax_region.shape

    # Collect all event-level information from each grid cell
    all_events = []

    for g in range(n_grid):
        res = detect_heatwaves_1d(
            years=years,
            max_temp=tasmax_region[g, :],
            min_temp=tasmin_region[g, :],
            max_threshold=max_threshold,
            min_threshold=min_threshold,
            min_consecutive_days=min_consecutive_days,
        )
        if len(res) == 0:
            continue

        res["Grid"] = g
        all_events.append(res)

    if not all_events:
        return pd.DataFrame(columns=[
            "Year", "Heatwave_Count", "Mean_Duration", "Mean_Intensity",
            "Region_Name", "Region_ID"
        ])

    all_events_df = pd.concat(all_events, ignore_index=True)

    # At this stage, each row corresponds to (Year, Grid) with metrics from detect_heatwaves_1d
    # First, aggregate by (Year, Grid): sum counts/duration/intensity within each grid cell
    grouped = all_events_df.groupby(["Year", "Grid"]).agg(
        Heatwave_Count=("Heatwave_Count", "sum"),
        # Note: using Mean_Duration here as a proxy for total duration.
        # store and sum the raw event durations separately.
        Total_Duration=("Mean_Duration", "sum"),
        Total_Intensity=("Mean_Intensity", "sum"),
    ).reset_index()

    # Then aggregate over all grid cells for each year
    agg = grouped.groupby("Year").agg(
        Heatwave_Count=("Heatwave_Count", "sum"),
        Duration=("Total_Duration", "sum"),
        Intensity=("Total_Intensity", "sum"),
    ).reset_index()

    total_grid_points = n_grid

    # Compute mean duration and intensity per heatwave (conditional on having events)
    agg["Mean_Duration"] = np.where(
        agg["Heatwave_Count"] > 0,
        agg["Duration"] / agg["Heatwave_Count"],
        0.0,
    )
    agg["Mean_Intensity"] = np.where(
        agg["Heatwave_Count"] > 0,
        agg["Intensity"] / agg["Heatwave_Count"],
        0.0,
    )

    # Heatwave_Count can be scaled to an average per grid cell
    agg["Heatwave_Count"] = agg["Heatwave_Count"] / total_grid_points

    agg["Region_Name"] = region_name
    agg["Region_ID"] = region_id

    return agg[["Year", "Heatwave_Count", "Mean_Duration", "Mean_Intensity",
                "Region_Name", "Region_ID"]]
def load_region_mask_indices(mask_path):
    """
    Build a dictionary of grid indices for each region based on
    'regionmask-region_osgb.nc'.

    The function extracts, for each region:
        - region_id
        - region_name
        - lists of latitude indices (LatID)
        - lists of longitude indices (LonID)

    Returns:
    dict
        Keys are region IDs (1-based). Each entry contains:
            {
                "Region_ID": int,
                "Region_Name": str,
                "LatID": array of y-indices,
                "LonID": array of x-indices,
            }
    """
    ds_mask = xr.open_dataset(mask_path)
    mask = ds_mask["region_mask"].values  # shape: (y, x, n_region)
    region_names = ds_mask["geo_region"].values.astype(str)

    # Apply name shortening
    region_names = np.char.replace(region_names, "North East", "NE")
    region_names = np.char.replace(region_names, "South East", "SE")
    region_names = np.char.replace(region_names, "North West", "NW")
    region_names = np.char.replace(region_names, "South West", "SW")
    region_names = np.char.replace(region_names, "Yorkshire and", "Yorks &")
    region_names = np.char.replace(region_names, "Northern", "N")
    region_names = np.char.replace(region_names, "North", "N")
    region_names = np.char.replace(region_names, "South", "S")
    region_names = np.char.replace(region_names, "East", "E")
    region_names = np.char.replace(region_names, "West", "W")

    regions_list = {}
    n_region = mask.shape[2]

    for r in range(n_region):
        ys, xs = np.where(mask[:, :, r] == 1)
        if ys.size == 0:
            continue
        regions_list[r + 1] = {
            "Region_ID": r + 1,
            "Region_Name": region_names[r],
            "LatID": ys,
            "LonID": xs,
        }

    ds_mask.close()
    return regions_list


def process_one_pair(tasmax_path, tasmin_path,
                     regions_list,
                     thresholds_df,
                     min_consecutive_days=3):
    """
    Process one pair of tasmax/tasmin NetCDF files and compute heatwave
    statistics for every region.

    This function:
        - extracts grid cells belonging to each region,
        - looks up the corresponding temperature thresholds,
        - applies the heatwave detection to each grid cell,
        - aggregates metrics across all grid cells in the region,
        - returns a combined DataFrame containing all regions.

    Parameters:
    tasmax_path : str, File path to tasmax NetCDF file.
    tasmin_path : str, File path to tasmin NetCDF file.
    regions_list : dict, Output from load_region_mask_indices(). Contains grid indices per region.
    thresholds_df : pandas.DataFrame
        Must contain: region, max_threshold, min_threshold.
    min_consecutive_days : int, Minimum number of consecutive days to qualify as a heatwave.

    Returns:
    pandas.DataFrame
        Rows contain regional heatwave statistics for all regions.
        Columns typically include:
            - Year
            - Heatwave_Count
            - Mean_Duration
            - Mean_Intensity
            - Region_Name
            - Region_ID
        If no results are found, returns an empty DataFrame.
    """
    ds_max = xr.open_dataset(tasmax_path)
    ds_min = xr.open_dataset(tasmin_path)

    tasmax = ds_max["tasmax"].values  # shape: (time, y, x)
    tasmin = ds_min["tasmin"].values
    years = ds_max["year"].values     # 'year' already provided in the files

    results = []

    for region_id, info in regions_list.items():
        region_name = info["Region_Name"]
        lat_ids = info["LatID"]
        lon_ids = info["LonID"]

        # Extract temperature time series for all grid cells in this region:
        # Each grid cell → a (time,) vector.
        # Stack into array of shape (n_grid, n_time).
        grid_series_max = []
        grid_series_min = []

        for lat, lon in zip(lat_ids, lon_ids):
            grid_series_max.append(tasmax[:, lat, lon])
            grid_series_min.append(tasmin[:, lat, lon])

        if len(grid_series_max) == 0:
            continue

        max_arr = np.stack(grid_series_max, axis=0)
        min_arr = np.stack(grid_series_min, axis=0)

        # Look up region-specific thresholds
        row = thresholds_df.loc[thresholds_df["region"] == region_name]
        if row.empty:
            print(f"Warning: no threshold found for region {region_name}")
            continue

        max_threshold = float(row["max_threshold"].iloc[0])
        min_threshold = float(row["min_threshold"].iloc[0])

        # Apply heatwave analysis for this region
        region_df = analyze_region_from_arrays(
            years=years,
            tasmax_region=max_arr,
            tasmin_region=min_arr,
            region_name=region_name,
            region_id=region_id,
            max_threshold=max_threshold,
            min_threshold=min_threshold,
            min_consecutive_days=min_consecutive_days,
        )

        results.append(region_df)

    if not results:
        return pd.DataFrame()

    ds_max.close()
    ds_min.close()

    return pd.concat(results, ignore_index=True)

