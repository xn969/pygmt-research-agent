"""Standalone recipe. Run: python figure.py path/to/spec.json OUTPUT_DIRECTORY"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pygmt
import xarray as xr
from pyproj import Geod


def read_grid(path, variable):
    with xr.open_dataset(path) as ds:
        grid = ds[variable].transpose("lat", "lon").load()
    for axis, limit in [("lat", 90), ("lon", 180)]:
        coords = grid[axis].values
        steps = np.diff(coords)
        if (not len(steps) or not np.isfinite(coords).all()
                or max(abs(coords)) > limit or min(steps) <= 0
                or not np.allclose(steps, steps[0], rtol=1e-4)):
            raise ValueError(f"{axis} must be regular, increasing geographic coordinates")
    if not np.isfinite(grid.values).any():
        raise ValueError(f"No finite samples in {path}")
    return grid.where(np.isfinite(grid))


def render(spec_path, out):
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    grids = []
    labels = []
    for item in spec["grids"]:
        grid = read_grid(spec_path.parent / item["file"], item["variable"])
        if not spec["unit"] or grid.attrs.get("units") != spec["unit"]:
            raise ValueError(f"Unit mismatch in {item['file']}; convert inputs first")
        grids.append(grid)
        labels.append(item["label"])
    if not 1 <= len(grids) <= 4 or len(set(labels)) != len(labels):
        raise ValueError("Supply one to four uniquely labelled grids")
    west = max(float(g.lon.min()) for g in grids)
    east = min(float(g.lon.max()) for g in grids)
    south = max(float(g.lat.min()) for g in grids)
    north = min(float(g.lat.max()) for g in grids)
    if west >= east or south >= north or max(abs(south), abs(north)) >= 85:
        raise ValueError("Regional Mercator maps require shared coverage within 85 degrees")
    region = [west, east, south, north]
    low, high = spec["scale"]
    if not np.isfinite([low, high]).all() or low >= high:
        raise ValueError("scale must be an increasing finite pair")
    fig = pygmt.Figure()
    pygmt.makecpt(cmap="viridis", series=[low, high])
    profile = None
    if "profile" in spec:
        line = spec["profile"]
        start, end = line["start"], line["end"]
        geod = Geod(ellps="WGS84")
        azimuth, _, length = geod.inv(*start, *end)
        step = line["step_km"]
        if not np.isfinite([step, length]).all() or min(step, length) <= 0:
            raise ValueError("Specify distinct endpoints and a positive step_km")
        intervals = int(np.ceil(length / (1000 * step)))
        if intervals > 100000:
            raise ValueError("Too many samples; increase step_km")
        distances = np.linspace(0, length, intervals + 1)
        ones = np.ones_like(distances)
        lon, lat, _ = geod.fwd(ones*start[0], ones*start[1], ones*azimuth, distances)
        if (lon.min() < west-1e-9 or lon.max() > east+1e-9
                or lat.min() < south-1e-9 or lat.max() > north+1e-9):
            raise ValueError("Profile extends outside common coverage")
        profile = pd.DataFrame(dict(lon=lon, lat=lat, distance_km=distances/1000))
        for grid, label in zip(grids, labels):
            if label in profile:
                raise ValueError(f"Reserved profile column: {label}")
            samples = grid.interp(lon=xr.DataArray(lon, dims="point"),
                                  lat=xr.DataArray(lat, dims="point"))
            if not np.isfinite(samples.values).any():
                raise ValueError(f"No valid profile values for {label}")
            profile[label] = samples.values
    for index, (grid, label) in enumerate(zip(grids, labels)):
        if index:
            fig.shift_origin(xshift="10c")
        fig.grdimage(grid=grid, region=region, projection="M8c", cmap=True,
                     frame=["af", f"WSne+t{label}"])
        if profile is not None:
            fig.plot(x=profile.lon, y=profile.lat, pen="1p,black")
    fig.colorbar(position="JMR+o0.6c/0c", frame=f"af+l{spec['unit']}")
    if profile is not None:
        fig.shift_origin(xshift=f"{-10*(len(grids)-1)}c", yshift="-7c")
        values = profile[labels].to_numpy()
        ymin, ymax = np.nanmin(values), np.nanmax(values)
        pad = max((ymax-ymin)*0.08, 1e-6)
        fig.basemap(region=[0, distances[-1]/1000, ymin-pad, ymax+pad],
                    projection="X12c/4c", frame=["xaf+lDistance (km)", f"yaf+l{spec['unit']}", "WSne"])
        for label, color in zip(labels, ["black", "red", "blue", "green"]):
            fig.plot(x=profile.distance_km, y=profile[label], pen=f"1p,{color}", label=label)
        fig.legend(position="JTR+jTR", box="+gwhite")
    out.mkdir(parents=True, exist_ok=False)
    fig.savefig(str(out / "figure.png"), dpi=180)
    if profile is not None:
        profile.to_csv(out / "profile.csv", index=False)
    with (out / "versions.txt").open("w", encoding="utf-8") as stream:
        pygmt.show_versions(file=stream)


if __name__ == "__main__":
    render(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
