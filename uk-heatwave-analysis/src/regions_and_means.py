import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


def load_region_mask(mask_path):
    """Load region mask NetCDF file and clean region names."""
    ds_mask = xr.open_dataset(mask_path)

    region_names = ds_mask["geo_region"].values.astype(str)
    mask = ds_mask["region_mask"].values  # shape: (y, x, n_region)
    osgb_x = ds_mask["projection_x_coordinate"].values
    osgb_y = ds_mask["projection_y_coordinate"].values

    # Rename the regions
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

    ds_mask.close()
    return mask, region_names, osgb_x, osgb_y


def plot_regions(mask, region_names, osgb_x, osgb_y, ukmap_path="ukmap2.dat"):
    """Plot region mask and overlay labels."""
    ukmap = pd.read_csv(ukmap_path, delim_whitespace=True)
    n_region = mask.shape[2]

    X, Y = np.meshgrid(osgb_x, osgb_y)

    plt.figure(figsize=(6, 8))
    plt.axis("equal")
    plt.axis("off")

    # Draw empty map outline first
    plt.plot(ukmap["East"], ukmap["North"], color="0.8", lw=2)
    
    for i in range(n_region):
        region_mask = mask[:, :, i]
        plt.contourf(X, Y, region_mask, levels=[0.5, 1.5], alpha=0.2)

    # Add region labels using centroid positions
    for i in range(n_region):
        region_mask = mask[:, :, i]
        ys, xs = np.where(region_mask == 1)
        if xs.size == 0:
            continue
        cx = osgb_x[xs].mean()
        cy = osgb_y[ys].mean()
        plt.text(cx, cy, region_names[i], ha="center", va="center", fontsize=8)

    plt.tight_layout()
    plt.savefig("RegionMap_python.png", dpi=150)
    # plt.show()


def compute_regional_means(tas_path, mask, region_names):
    """
    Compute regional average temperatures
    For each region and each time step, compute the spatial average.
    """
    ds = xr.open_dataset(tas_path)
    tas = ds["tas"]  # shape: (time, y, x)
    time = ds["time"]

    # Broadcast mask from (y, x, region) to (time, y, x, region)
    mask_da = xr.DataArray(mask, dims=("y", "x", "region"))
    tas_expanded = tas.expand_dims({"region": mask_da.region}).transpose(
        "time", "y", "x", "region"
    )

    mask_bool = mask_da > 0

    # Apply mask per region: only keep grid cells inside region
    # Then average over spatial dimensions (y, x)
    regional_means = tas_expanded.where(mask_bool).mean(
        dim=("projection_y_coordinate", "projection_x_coordinate"),
        skipna=True
    )

    dates = xr.decode_cf(ds).time.to_pandas()
    df = pd.DataFrame(index=dates)

    # Build DataFrame with one column per region (16 Regions in total)
    for i, name in enumerate(region_names):
        df[name.replace(" ", "_")] = regional_means.isel(region=i).values

    ds.close()
    return df.reset_index().rename(columns={"index": "Date"})


def main():
    """
    Example usage of the region-masking and temperature-aggregation functions.

    This script:
      1. Loads the OSGB region mask
      2. Plots the UK regional boundaries
      3. Computes regional mean temperatures from HadUK-Grid
      4. Saves the results as a CSV file

    Note:
        This function is only provided as an example workflow.
        Users may call the functions individually in other scripts.
    """
    tas_path = "tas_hadukgrid_uk_12km_mon_188401-202112.nc"
    mask_path = "regionmask-region_osgb.nc"

    mask, region_names, osgb_x, osgb_y = load_region_mask(mask_path)
    plot_regions(mask, region_names, osgb_x, osgb_y)

    regional_df = compute_regional_means(tas_path, mask, region_names)
    regional_df.to_csv("tas_regional_python.csv", index=False)



if __name__ == "__main__":
    main()
