import unittest

from pygmt_research.figure import annotation_sides


class TestAxes(unittest.TestCase):
    def test_maps_over_profile(self):
        panels = [(0, 0, 1, 1, "lon", "lat"), (1, 0, 2, 1, "lon", "lat"),
                  (0, 1, 2, 2, "distance", "velocity")]
        self.assertEqual(annotation_sides(panels), ["WseN", "wsEN", "WSEn"])

    def test_enclosed_scales(self):
        panels = [(x, y, x+1, y+1, "lon", "lat") for y in range(3) for x in range(3)]
        self.assertEqual(annotation_sides(panels)[4], "wsen")
        panels[4] = (1, 1, 2, 2, "time", "depth")
        self.assertEqual(annotation_sides(panels)[4], "WseN")
