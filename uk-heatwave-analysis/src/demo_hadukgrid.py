import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from netCDF4 import Dataset


def get_nc_var_info(path):
    ds = Dataset(path)
    info = []
    for name, var in ds.variables.items():
        units = getattr(var, "units", None)
        info.append({"Name": name, "Units": units})
    ds.close()
    return pd.DataFrame(info)


def main():
    # ---- 1. Read UK coastline (ukmap2.dat) ----
    # Assuming it is whitespace-delimited 
    ukmap = pd.read_csv("ukmap2.dat", delim_whitespace=True)

    # ---- 2. Open the HadUK-Grid file ----
    nc_path = "tas_hadukgrid_uk_12km_mon_188401-202112.nc"

    print(get_nc_var_info(nc_path))

    ds = xr.open_dataset(nc_path)

    # Variable name may vary depending on file; “tas” matches the R script
    tas = ds["tas"]   # Usually shaped (time, y, x) or (time, projection_y, projection_x)

    print("Dimensions of temperature data:", tas.shape)

    # Coordinates
    time = ds["time"]
    osgb_x = ds["projection_x_coordinate"]
    osgb_y = ds["projection_y_coordinate"]

    print("Length of time variable:", time.size)
    print("Length of OSGB.x:", osgb_x.size)
    print("Length of OSGB.y:", osgb_y.size)

    # ---- 3. Compute spatial and temporal means ----
    # R: MeanTemp <- rowMeans(HadUKTemp, dims=2) → mean over time → (y, x)
    mean_temp = tas.mean(dim="time")  # DataArray: (y, x)

    # R: UKMean <- colMeans(HadUKTemp, dims=2)  → mean over all grid cells → (time)
    uk_mean = tas.mean(dim=("projection_y_coordinate", "projection_x_coordinate"))

    # ---- 4. Plot spatial distribution of UK mean temperature ----
    plt.figure(figsize=(6, 8))

    # Ensure meshgrid order matches mean_temp’s dimensions
    X, Y = np.meshgrid(osgb_x.values, osgb_y.values)
    plt.pcolormesh(X, Y, mean_temp.values, shading="auto")
    plt.colorbar(label="Mean temperature (°C)")
    plt.title("Mean temperature (°C)")
    plt.axis("equal")

    # Overlay coastline
    plt.plot(ukmap["East"], ukmap["North"], linewidth=2, color="k")

    plt.tight_layout()
    plt.savefig("UKMeanTemp_python.png", dpi=150)
    # plt.show()

    # ---- 5. Process time variable: days since 1884-01-16 ----
    # If xarray already decodes CF time, you can use time.dt.year.
    # If not, decode manually based on units.
    # Here we assume time is an offset in days.
    time0 = np.datetime64("1884-01-16")
    dates = time0 + time.values.astype("timedelta64[D]")
    years = pd.DatetimeIndex(dates).year

    # Build the time series for mean UK temperature
    uk_mean_series = pd.Series(uk_mean.values, index=pd.DatetimeIndex(dates))

    plt.figure(figsize=(10, 4))
    plt.plot(
        uk_mean_series.index.year + (uk_mean_series.index.month - 0.5) / 12.0,
        uk_mean_series.values,
        lw=0.8
    )
    plt.xlabel("Year")
    plt.ylabel("°C")
    plt.title("Mean monthly UK temperature")

    overall_mean = np.nanmean(uk_mean_series.values)
    plt.axhline(overall_mean, ls="--", lw=2)

    plt.tight_layout()
    plt.savefig("UKMeanTemp_timeseries_python.png", dpi=150)
    # plt.show()

    ds.close()


if __name__ == "__main__":
    main()
