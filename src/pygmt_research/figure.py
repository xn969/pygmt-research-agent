"""Standalone recipe. Run: python figure.py path/to/spec.json OUTPUT_DIRECTORY"""
import json
import os
import re
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pygmt
import xarray as xr
from pyproj import Geod
from PIL import Image
from fontTools.agl import UV2AGL
from fontTools.pens.basePen import BasePen
from fontTools.ttLib import TTFont


def annotation_sides(panels):
    """Logical rectangles are (left, top, right, bottom, x_axis, y_axis)."""
    marked = [set() for _ in panels]
    for axis, before, after in [(0, "N", "S"), (1, "W", "E")]:
        cross = 1 - axis
        connected = [set() for _ in panels]
        for i, panel in enumerate(panels):
            for side, edge, opposite in [(before, cross, cross+2), (after, cross+2, cross)]:
                neighbors = [j for j, other in enumerate(panels)
                             if other[opposite] == panel[edge] and
                             max(panel[axis], other[axis]) < min(panel[axis+2], other[axis+2])]
                if not neighbors:
                    marked[i].add(side)
                for j in neighbors:
                    other = panels[j]
                    if (panel[axis], panel[axis+2], panel[4+axis]) == (
                            other[axis], other[axis+2], other[4+axis]):
                        connected[i].add(j)
        # An enclosed group with its own scale gets one consistent annotation.
        unseen = set(range(len(panels)))
        while unseen:
            stack = [unseen.pop()]
            group = []
            while stack:
                i = stack.pop()
                group.append(i)
                pending = connected[i] & unseen
                stack.extend(pending)
                unseen.difference_update(pending)
            if not any(marked[i] & {before, after} for i in group):
                marked[min(group, key=lambda i: panels[i][cross])].add(before)
    return ["".join(s if s in sides else s.lower() for s in "WSEN") for sides in marked]


def check_layout(path, boxes, gap, dpi):
    """Check expected panel corners against horizontal frame rules in the PNG."""
    with Image.open(path) as image:
        ink = np.asarray(image.convert("RGB")).max(axis=2) < 65
    px = dpi / 2.54
    tolerance = 7
    rules = []
    min_width = min(box[2]-box[0] for box in boxes)*px-tolerance
    for y, row in enumerate(ink):
        edges = np.flatnonzero(np.diff(np.r_[False, row, False])).reshape(-1, 2)
        rules.extend((int(left), y, int(right)) for left, right in edges
                     if right-left >= min_width)
    first_width = (boxes[0][2]-boxes[0][0])*px
    origins = [rule for rule in rules if abs(rule[2]-rule[0]-first_width) <= tolerance]
    for x0, y0, _ in origins:
        measured = []
        for left, top, right, bottom in boxes:
            def locate(y):
                target = np.array([x0+left*px, y0+y*px, x0+right*px])
                nearest = min(rules, key=lambda r: max(abs(r-target)))
                return nearest if max(abs(nearest-target)) <= tolerance else None
            upper, lower = locate(top), locate(bottom)
            if upper is None or lower is None:
                break
            measured.append([upper[0], upper[1], upper[2], lower[1]])
        if len(measured) == len(boxes):
            return {"frame_bounds_px": measured, "gutter_cm": gap, "tolerance_px": tolerance}
    raise ValueError("Layout check: expected panel frames could not be verified")


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
             text=f"({tag})", font="12p,Times-Bold,black",
             fill="white@30", clearance="0.05c/0.04c")


class OutlinePen(BasePen):
    """Write a TrueType contour as a PostScript path, including cubic conversion."""
    def __init__(self, glyphs):
        super().__init__(glyphs)
        self.commands = []

    def command(self, name, *points):
        self.commands.append(" ".join(f"{v:.3f}" for point in points for v in point) + " " + name)

    def _moveTo(self, point):
        self.command("moveto", point)

    def _lineTo(self, point):
        self.command("lineto", point)

    def _curveToOne(self, *points):
        self.command("curveto", *points)

    def _closePath(self):
        self.commands.append("closepath")


def outline_font(file, name):
    # Keep the artwork portable across Ghostscript versions. The temporary EPS
    # contains outlines for the Western PostScript glyph names used by GMT.
    with TTFont(file) as font:
        glyphs = font.getGlyphSet()
        head = font["head"]
        scale = 1 / head.unitsPerEm
        result = [f"/{name} 12 dict begin /FontType 3 def /FontName /{name} def",
                  f"/FontMatrix [{scale} 0 0 {scale} 0 0] def",
                  f"/FontBBox [{head.xMin} {head.yMin} {head.xMax} {head.yMax}] def",
                  "/Encoding StandardEncoding def /Glyphs 800 dict dup begin",
                  "/.notdef {0 0 setcharwidth} bind def"]
        for code, glyph in font.getBestCmap().items():
            if code not in UV2AGL:
                continue
            pen = OutlinePen(glyphs)
            glyphs[glyph].draw(pen)
            result.append(f"/{UV2AGL[code]} {{{glyphs[glyph].width} 0 setcharwidth newpath "
                          + " ".join(pen.commands) + " fill} bind def")
        result.extend(["end def /BuildGlyph {exch /Glyphs get exch",
                       "2 copy known not {pop /.notdef} if get exec} bind def",
                       "/BuildChar {1 index /Encoding get exch get",
                       "1 index /BuildGlyph get exec} bind def currentdict end definefont pop"])
        return "\n".join(result).encode("ascii")


def export_png(fig, path):
    """Map GMT's Times font slots to the user's installed TrueType family."""
    roots = [Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts",
             Path("/Library/Fonts/Microsoft"), Path("/System/Library/Fonts/Supplemental"),
             Path("/usr/share/fonts/truetype/msttcorefonts")]
    if os.getenv("PYGMT_RESEARCH_FONT_DIR"):
        roots.insert(0, Path(os.environ["PYGMT_RESEARCH_FONT_DIR"]))
    for root in roots:
        for regular, bold in [("times.ttf", "timesbd.ttf"),
                              ("Times New Roman.ttf", "Times New Roman Bold.ttf"),
                              ("Times_New_Roman.ttf", "Times_New_Roman_Bold.ttf")]:
            files = [root / regular, root / bold]
            if not all(f.is_file() for f in files):
                continue
            with tempfile.TemporaryDirectory() as tmp:
                eps = Path(tmp) / "figure.eps"
                fig.savefig(str(eps))
                ps = eps.read_bytes().replace(b"/Times-Roman ", b"/TimesNewRomanPSMT ")
                ps = ps.replace(b"/Times-Bold ", b"/TimesNewRomanPS-BoldMT ")
                fonts = b"\n".join(outline_font(f, n) for f, n in
                                   zip(files, ["TimesNewRomanPSMT", "TimesNewRomanPS-BoldMT"]))
                eps.write_bytes(ps.replace(b"%%BeginProlog", b"%%BeginProlog\n" + fonts, 1))
                subprocess.run(["gmt", "psconvert", str(eps), "-Tg", "-A", "-E180",
                                "-Qt4", "-Qg4", f"-F{path.with_suffix('')}"], check=True,
                               cwd=tmp, env={**os.environ, "GMT_SESSION_NAME": "eps-export"})
            return "Times New Roman"
    warnings.warn("Times New Roman unavailable; using the Ghostscript Times substitute.")
    fig.savefig(str(path), dpi=180)
    return "Times substitute"


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
    pygmt.config(MAP_FRAME_TYPE="plain", FONT="11p,Times-Roman",
                 FONT_ANNOT_PRIMARY="9p,Times-Roman", FORMAT_GEO_MAP="ddd.xx")
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
    columns = spec.get("columns", len(grids))
    if type(columns) is not int or not 1 <= columns <= len(grids):
        raise ValueError("columns must be between 1 and the number of grids")
    rows = (len(grids) + columns - 1) // columns
    span = width + stride * (columns - 1)
    bounds = "/".join(str(v) for v in region)
    map_height = float(subprocess.check_output(
        ["gmt", "mapproject", f"-R{bounds}", f"-JM{width}c", "-Wh", "--PROJ_LENGTH_UNIT=cm"],
        text=True))
    panels, boxes = [], []
    for i in range(len(grids)):
        row, col = divmod(i, columns)
        panels.append((col, row, col+1, row+1, ("longitude", west, east),
                       ("mercator latitude", south, north)))
        left, top = col*stride, row*(map_height+gap)
        boxes.append((left, top, left+width, top+map_height))
    if profile is not None:
        panels.append((0, rows, columns, rows+1, "distance_km", ("value", spec["unit"])))
        top = rows*(map_height+gap)
        boxes.append((0, top, span, top+profile_height))
    frames = annotation_sides(panels)
    for index, (grid, label) in enumerate(zip(grids, labels)):
        if index:
            fig.shift_origin(xshift=f"{boxes[index][0]-boxes[index-1][0]}c",
                             yshift=f"{boxes[index-1][3]-boxes[index][3]}c")
        sides = frames[index]
        fig.grdimage(grid=grid, region=region, projection=f"M{width}c", cmap=True,
                     frame=["af", f"{sides}+t{label}" if spec.get("panel_titles", False) else sides])
        if profile is not None:
            fig.plot(x=profile.lon, y=profile.lat, pen="1p,black")
        mark_panel(fig, tags[index])
        if index == columns-1:
            # Plot coordinates are relative to the lower left of the top-right map.
            group_height = rows * map_height + (rows - 1) * gap
            center_y = map_height - group_height / 2
            fig.colorbar(position=f"x{width + 1.3}c/{center_y}c+jML"
                                  f"+w{0.8 * map_height}c/0.25c+v",
                         frame=f"af+l{spec['unit']}")
    if profile is not None:
        fig.shift_origin(xshift=f"{-boxes[-2][0]}c",
                         yshift=f"{-(profile_height + gap)}c")
        values = profile[labels].to_numpy()
        ymin, ymax = np.nanmin(values), np.nanmax(values)
        pad = max((ymax-ymin)*0.08, 1e-6)
        fig.basemap(region=[0, distances[-1]/1000, ymin-pad, ymax+pad],
                    projection=f"X{span}c/{profile_height}c",
                    frame=["xaf+lDistance (km)", f"yaf+l{spec['unit']}", frames[-1]])
        for label, color in zip(labels, ["black", "red", "blue", "green"]):
            fig.plot(x=profile.distance_km, y=profile[label], pen=f"1p,{color}", label=label)
        fig.legend(position="JTR+jTR", box="+gwhite")
        mark_panel(fig, tags[-1])
    out.mkdir(parents=True, exist_ok=False)
    family = export_png(fig, out / "figure.png")
    report = check_layout(out / "figure.png", boxes, gap, 180)
    report["annotated_sides"] = ["".join(s for s in frame if s.isupper()) for frame in frames]
    report["tags"] = tags
    report["font"] = family
    (out / "layout.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    if profile is not None:
        profile.to_csv(out / "profile.csv", index=False)
    with (out / "versions.txt").open("w", encoding="utf-8") as stream:
        pygmt.show_versions(file=stream)


if __name__ == "__main__":
    render(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
