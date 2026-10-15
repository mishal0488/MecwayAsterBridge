"""Beam section forces retain element identity, local axes and time steps."""
import os
import sys
import unittest
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'Scripts'))
import aster_results_export as export


class Fields:
    n = 4
    rows = None
    results = {'res': object()}

    def __init__(self):
        self.layouts = {('res', 'SIEF_ELNO'):
                        (np.array([0, 2, 4]), np.array([2, 2, 3]),
                         np.array([1, 1, 1]))}
        self.values = np.arange(42, dtype=float).reshape(7, 6)
        self.mask = np.ones((7, 6), dtype=bool)

    def get(self, rkey, field, on_cells):
        if field != 'SIEF_ELNO':
            return None
        return (['N', 'VY', 'VZ', 'MT', 'MFY', 'MFZ'],
                self.values, self.mask, None)


class BeamResultantsTests(unittest.TestCase):
    def setUp(self):
        self.fields = Fields()
        # A shared node, reversed solver connectivity, and a shell with the
        # same component names: none may mix the local beam quantities.
        self.conn = [[1, 0], [1, 2], [0, 2, 3]]
        self.model = {'elements': {10: ('line2', [1, 2]),
                                   20: ('line2', [2, 3]),
                                   30: ('tri3', [1, 3, 4])}}

    def plan(self, groups=None):
        return export._plan_variables(self.fields, self.model, self.conn, groups)

    def test_native_names_axes_and_discontinuous_junction(self):
        variables = self.plan()
        self.assertEqual([v.name for v in variables],
                         ['tensileforce', 'sheary', 'shearz', 'torsion',
                          'bendingmomv', 'bendingmomw'])
        for col, var in enumerate(variables):
            self.assertTrue(var.element)
            np.testing.assert_array_equal(var.ids, [10, 20])
            np.testing.assert_array_equal(var.getter(self.fields),
                                          self.fields.values[[1, 0, 2, 3], col])

    def test_later_time_values_are_read_again(self):
        variables = self.plan()
        self.fields.values = self.fields.values * 0.25
        for col, var in enumerate(variables):
            np.testing.assert_array_equal(var.getter(self.fields),
                                          self.fields.values[[1, 0, 2, 3], col])

    def test_missing_component_does_not_export_invalid_element(self):
        self.fields.mask[0, 3] = False
        variables = {v.name: v for v in self.plan()}
        np.testing.assert_array_equal(variables['torsion'].ids, [20])
        np.testing.assert_array_equal(variables['sheary'].ids, [10, 20])

    def test_result_filter(self):
        self.assertEqual(self.plan(['displacements']), [])
        self.assertEqual(len(self.plan(['stresses'])), 6)


if __name__ == '__main__':
    unittest.main()
