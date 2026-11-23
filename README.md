# PyGMT Research Agent

Make regional raster comparisons and WGS84 geodesic profiles from your own NetCDF files. Includes a command-line tool and a portable Agent Skill. The examples generate synthetic data locally.

## Try it

From the extracted source folder, with Miniforge installed:

```sh
conda env create -f environment.yml
conda activate pygmt-research
python -m pip install --no-deps .
pygmt-research doctor
pygmt-research demo comparison runs/comparison
pygmt-research demo profile runs/profile
```

Open `runs/profile/output/figure.png`. The adjacent profile.csv contains coordinates, distance in km and interpolated values. versions.txt identifies the runtime. Edit the generated spec.json and run:

```sh
pygmt-research render runs/profile/spec.json runs/revised
python runs/profile/output/figure.py runs/profile/spec.json runs/independent
```

Use a new output directory for each run. The standalone script needs the same Python/GMT environment and original input grids. For desktop agents invoke the CLI through `conda run -n pygmt-research pygmt-research ...` if the shell isn't activated. Copy the skills/pygmt-research folder to a skill discovery directory supported by your host. Try asking “Compare my two temperature grids using one scale and sample both along this transect.”

## Input contract

The generated spec is a working example. Set grids to one to four objects with file (relative to the spec or absolute), variable and a unique short label. Each variable must be two-dimensional with regular increasing lat and lon axes in WGS84 degrees, within ±180 longitude. The units attribute must match spec.unit. Select a shared increasing scale pair. Panels use the intersection of grid extents; use regional extents below 85 degrees latitude and within one longitude interval. Convert projected grids and harmonize physical units beforehand.

Omit profile for maps only. For profile set start and end as [longitude, latitude] and step_km as a maximum sampling interval. Points follow the WGS84 geodesic, including both endpoints. Each grid is bilinearly interpolated on its native axes, and missing samples remain missing. No subtraction, smoothing or reference-date adjustment is implicit.

## Panel layout

Set `"columns": 2` for two map columns; by default all maps occupy one row.
Coordinates appear on exposed panel edges. Aligned neighbors sharing a
coordinate scale use outer labels. When an enclosed group has a distinct
scale, it uses the top for x or the left for y. A full-width profile uses
distance and data values independently of map longitude and latitude.
The layout report records annotated sides (W, S, E, N) in panel order.

Maps share a frame size and a compact 0.85 cm gutter. The optional profile
spans the full map row with the same frame-to-frame gutter. A shared colorbar
sits outside the last map. Panel order is left to right, then the profile.
In spec.json, set `"panel_labels": "lower"` for (a), (b), or use `"numeric"`
or a list such as `["B", "C", "D"]`. The default is upper case. Tags have
black lettering and a white background with 30% transparency.

Map titles are hidden by default. Set `"panel_titles": true` to show the
grid labels above the maps. Labels still identify profile series and CSV columns.
Figures use locally installed Times New Roman (regular and bold). For a
nonstandard font location set PYGMT_RESEARCH_FONT_DIR to the font directory.
When unavailable, rendering warns and records the Times substitute in layout.json.

The exported PNG is checked for aligned frame rules and equal gutters.
See layout.json for measured panel bounds in pixels.

## License

MIT for this project's code. Dependencies retain their own licenses.
