import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from geodata.model.wind import wind_speed_merra2


def main():
    """Main function."""
    ds = xr.open_dataset(
        "/home/xiqiang/cli/MERRA2_300.tavg1_2d_slv_flx_Nx.20100101.nc4"
    )

    ds = wind_speed_merra2(ds, compute_lml=False)

    alpha = ds["coeffs"][..., 0].values[..., np.newaxis]
    beta = ds["coeffs"][..., 1].values[..., np.newaxis]

    heights = np.linspace(2, 140)
    log_heights = np.log(heights - ds["disph"].values[..., np.newaxis])
    speeds = alpha * log_heights + beta

    c2m = (ds["u2m"].values[2, 30, 30] ** 2 + ds["v2m"].values[2, 30, 30] ** 2) ** 0.5
    c10m = (
        ds["u10m"].values[2, 30, 30] ** 2 + ds["v10m"].values[2, 30, 30] ** 2
    ) ** 0.5
    c50m = (
        ds["u50m"].values[2, 30, 30] ** 2 + ds["v50m"].values[2, 30, 30] ** 2
    ) ** 0.5
    clml = (
        ds["ulml"].values[2, 30, 30] ** 2 + ds["vlml"].values[2, 30, 30] ** 2
    ) ** 0.5
    hlml = ds["hlml"].values[2, 30, 30]

    print([c2m, c10m, c50m, clml, ds["disph"].values[2, 30, 30]])
    print(alpha[2, 30, 30], beta[2, 30, 30])
    print(speeds[2, 30, 30])

    plt.plot(speeds[2, 30, 30], heights)
    plt.plot([c2m, c10m, c50m, clml], [2, 10, 50, hlml], "o")
    plt.ylabel("Height (m)")
    plt.xlabel("Wind Speed (m/s)")
    plt.title("Wind Speed vs Height at [2, 30, 30]")
    plt.tight_layout()
    plt.savefig("speed.png")


if __name__ == "__main__":
    main()
