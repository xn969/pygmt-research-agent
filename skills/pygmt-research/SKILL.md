---
name: pygmt-research
description: Create publication figures comparing regional geographic grids, with optional geodesic profiles, using PyGMT. Use for scientific map panels, shared scales and along-line comparisons.
---

Use the user's configured plotting environment. Run `pygmt-research doctor` there when first setting up. For a synthetic reference, run `pygmt-research demo profile NEW_DIRECTORY` (or comparison) and inspect the JSON and image.

Inspect the user's NetCDF variable names, coordinate axes, coverage, units and missing values before assembling a spec. Use explicit file-to-panel mappings, consistent physical units and one agreed scale. Resolve ambiguous units, reference epochs or subtraction direction with the researcher. This recipe plots supplied quantities directly.

Run `pygmt-research render SPEC.json NEW_OUTPUT_DIR`. Read the status, inspect the PNG and check sampled distances and finite values. Report the figure, editable script, spec and profile table if present. Warn about saturation if values exceed the selected scale. Keep shared scientific parameters fixed across a batch.

For an unfamiliar layout, start with a copy of the editable figure.py and adapt native PyGMT code in the user's project. Preserve the original inputs and the validated recipe. Record any changes to sampling or scientific assumptions.

Use panel_labels in the spec for user-requested lettering (upper, lower, numeric or a list). Read layout.json after rendering and inspect the preview for long titles, clipped annotations and legend overlap. Keep the full-width profile aligned with the outer map frames.

Use Times New Roman and omit individual map titles by default. Set panel_titles
to true only when requested. Mention a reported font substitution when delivering a figure.
