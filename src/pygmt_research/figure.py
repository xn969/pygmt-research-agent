"""Standalone recipe. Run: python figure.py path/to/spec.json OUTPUT_DIRECTORY"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pygmt
import xarray as xr
from pyproj import Geod
from PIL import Image


def check_layout(path, columns, profile, width, gap, height, dpi):
    """Locate frame rules in the exported PNG and verify their geometry."""
    with Image.open(path) as image:
        ink = np.asarray(image.convert("RGB")).max(axis=2) < 65
    px = dpi / 2.54
    panels = []
    tolerance = 7  # Includes both rasterized stroke edges at 180 dpi.

    def matching_rules(count, frame_width):
        matches = []
        for y, row in enumerate(ink):
            edges = np.flatnonzero(np.diff(np.r_[False, row, False]))
            runs = edges.reshape(-1, 2)
            runs = runs[np.abs(runs[:, 1] - runs[:, 0] - frame_width * px) <= tolerance]
            if len(runs) == count:
                matches.append((y, runs))
        return matches

    maps = matching_rules(columns, width)
    if not maps:
        raise ValueError("Layout check: could not locate aligned map frames")
    top, bounds = maps[0]
    # A single-column profile shares the map width, so find its rules by height.
    bottom = maps[-1][0]
    if profile:
        span = columns * width + (columns - 1) * gap
        rules = matching_rules(1, span)
        if not rules:
            raise ValueError("Layout check: missing profile frame")
        base = rules[-1][0]
        target = base - height * px
        row, profile_bounds = min(rules, key=lambda entry: abs(entry[0] - target))
        candidates = [y for y, _ in maps if y < row - tolerance]
        if not candidates:
            raise ValueError("Layout check: missing lower map frame")
        bottom = max(candidates)
        errors = [row-target, row-bottom-gap*px, profile_bounds[0, 0]-bounds[0, 0],
                  profile_bounds[0, 1]-bounds[-1, 1]]
        if max(map(abs, errors)) > tolerance:
            raise ValueError("Layout check: profile alignment or vertical gutter differs")
        panels.append([int(bounds[0, 0]), row, int(bounds[-1, 1]), base])
    if bottom <= top or np.any(abs(bounds[1:, 0]-bounds[:-1, 1]-gap*px) > tolerance):
        raise ValueError("Layout check: map frames or gutters differ")
    panels[:0] = [[int(left), top, int(right), bottom] for left, right in bounds]
    return {"frame_bounds_px": panels, "gutter_cm": gap, "tolerance_px": tolerance}


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


def panel_tags(style, count):
    if isinstance(style, list):
        if (len(style) != count or any(not isinstance(s, str) or not
                re.fullmatch(r"[A-Za-z0-9]{1,8}", s) for s in style)
                or len(set(style)) != count):
            raise ValueError("Supply one unique 1-8 character alphanumeric tag per panel")
        return style
    if style not in ("upper", "lower", "numeric"):
        raise ValueError("panel_labels must be upper, lower, numeric or a list")
    if style == "numeric":
        return [str(i + 1) for i in range(count)]
    start = ord("A" if style == "upper" else "a")
    return [chr(start + i) for i in range(count)]


def mark_panel(fig, tag):
    fig.text(position="TL", justify="TL", offset="0.16c/-0.16c",
             text=f"({tag})", font="12p,Helvetica-Bold,black",
             fill="white@30", clearance="0.05c/0.04c")


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
    pygmt.config(MAP_FRAME_TYPE="plain", FONT_ANNOT_PRIMARY="9p", FORMAT_GEO_MAP="ddd.xx")
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
    tags = panel_tags(spec.get("panel_labels", "upper"), len(grids) + int(profile is not None))
    width, gap, profile_height = 8, 0.85, 4
    stride = width + gap
    span = width + stride * (len(grids) - 1)
    for index, (grid, label) in enumerate(zip(grids, labels)):
        if index:
            fig.shift_origin(xshift=f"{stride}c")
        sides = "WSne" if index == 0 else "wSne"
        fig.grdimage(grid=grid, region=region, projection=f"M{width}c", cmap=True,
                     frame=["af", f"{sides}+t{label}"])
        if profile is not None:
            fig.plot(x=profile.lon, y=profile.lat, pen="1p,black")
        mark_panel(fig, tags[index])
    fig.colorbar(position="JMR+o0.6c/0c", frame=f"af+l{spec['unit']}")
    if profile is not None:
        fig.shift_origin(xshift=f"{-stride*(len(grids)-1)}c",
                         yshift=f"{-(profile_height + gap)}c")
        values = profile[labels].to_numpy()
        ymin, ymax = np.nanmin(values), np.nanmax(values)
        pad = max((ymax-ymin)*0.08, 1e-6)
        fig.basemap(region=[0, distances[-1]/1000, ymin-pad, ymax+pad],
                    projection=f"X{span}c/{profile_height}c",
                    frame=["xaf+lDistance (km)", f"yaf+l{spec['unit']}", "WSne"])
        for label, color in zip(labels, ["black", "red", "blue", "green"]):
            fig.plot(x=profile.distance_km, y=profile[label], pen=f"1p,{color}", label=label)
        fig.legend(position="JTR+jTR", box="+gwhite")
        mark_panel(fig, tags[-1])
    out.mkdir(parents=True, exist_ok=False)
    fig.savefig(str(out / "figure.png"), dpi=180)
    report = check_layout(out / "figure.png", len(grids), profile is not None,
                          width, gap, profile_height, 180)
    report["tags"] = tags
    (out / "layout.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    if profile is not None:
        profile.to_csv(out / "profile.csv", index=False)
    with (out / "versions.txt").open("w", encoding="utf-8") as stream:
        pygmt.show_versions(file=stream)


if __name__ == "__main__":
    render(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
