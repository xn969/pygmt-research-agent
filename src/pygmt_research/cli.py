"""Generate a trial or run a user specification."""
import argparse
import json
from pathlib import Path
import shutil
import sys
import tempfile


def demo(folder, profile):
    import numpy as np
    import xarray as xr

    folder.mkdir(parents=True, exist_ok=False)
    lon = np.linspace(10, 11, 81)
    lat = np.linspace(44, 44.7, 61)
    x, y = np.meshgrid(lon-10.5, lat-44.35)
    spec = {"unit": "degC", "scale": [8, 22], "grids": []}
    for index, label in enumerate(["Period A", "Period B"]):
        values = 10 + 10*np.exp(-18*((x-index*0.12)**2+y**2))
        grid = xr.DataArray(values, coords={"lat": lat, "lon": lon},
                            dims=("lat", "lon"), name="temperature", attrs={"units": "degC"})
        filename = f"period-{index+1}.nc"
        grid.to_netcdf(folder / filename)
        spec["grids"].append({"file": filename, "variable": grid.name, "label": label})
    if profile:
        spec["profile"] = {"start": [10.1, 44.15], "end": [10.9, 44.55], "step_km": 1}
    path = folder / "spec.json"
    path.write_text(json.dumps(spec, indent=2), encoding="utf-8")
    return path


def main():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="action", required=True)
    commands.add_parser("doctor")
    example = commands.add_parser("demo")
    example.add_argument("kind", choices=["comparison", "profile"])
    example.add_argument("folder", type=Path)
    run = commands.add_parser("render")
    run.add_argument("spec", type=Path)
    run.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        from .figure import render
        if args.action == "doctor":
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                render(demo(root / "data", False), root / "output")
                print((root / "output/versions.txt").read_text(encoding="utf-8"))
            return
        if args.action == "demo":
            spec = demo(args.folder, args.kind == "profile")
            output = args.folder / "output"
        else:
            spec, output = args.spec, args.output
        render(spec.resolve(), output.resolve())
        shutil.copyfile(Path(__file__).with_name("figure.py"), output / "figure.py")
        print(output.resolve() / "figure.png")
    except (OSError, ValueError, KeyError, ImportError, RuntimeError) as error:
        parser.exit(1, f"{error}\n")


if __name__ == "__main__":
    main()
