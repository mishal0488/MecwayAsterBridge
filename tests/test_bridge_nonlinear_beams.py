# -*- coding: utf-8 -*-
"""Offline end-to-end translation tests; no Code_Aster or Mecway needed."""
import copy
import io
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "Scripts"))
sys.path.insert(0, os.path.join(ROOT, "TestFiles", "Verificatie"))
import bridge_utils as utils
import bridge_from_liml as parser
import bridge_mesh as mesh
import bridge_to_aster as writer
from verif_lib import Model, linemesh


def beam_model(section=None, plastic=False, distributed=False):
    m = Model("L30")
    m.elset("Beam", "Steel")
    eids, nodes = linemesh(m, "Beam", 8, (0, 0, 0), (1, 0, 0))
    m.mat("Steel", section or
          '<geometric type="RectBar" recta="0.02 m" rectb="0.03 m" />',
          '<mechanical type="Isotropic" youngsmodulus="210000 MPa" '
          'poissonratio="0.3" />',
          '<density density="7850 kg/m3" />',
          '<plastic type="BilinearIsotropic" yieldstrength="250 MPa" '
          'tangentmodulus="1000 MPa" />' if plastic else "")
    m.node_selection("Root", [nodes[0]])
    m.node_selection("Tip", [nodes[-1]])
    m.add('<fixed selection="Root" />')
    m.add('<force selection="Tip"><x>0 N</x><y time="0 s">0 N</y>'
          '<y time="1 s">100 N</y><z>0 N</z></force>')
    if distributed:
        m.element_selection("Span", eids)
        m.add('<linepressure selection="Span"><x>0 N/m</x>'
              '<y>10 N/m</y><z>0 N/m</z></linepressure>')
    return m


class BeamTranslationTests(unittest.TestCase):
    def setUp(self):
        self.config = copy.deepcopy(utils.CONFIG_BRIDGE)
        utils.CONFIG_BRIDGE.update(active_analysis="nonlin_static_3d",
                                   elem_beam="POU_D_T_GD", med_export=True)
        self.logs = []
        self.loggers = [(mod, mod.log) for mod in (utils, parser, mesh, writer)]
        for mod, _ in self.loggers:
            mod.log = self.logs.append
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()
        utils.CONFIG_BRIDGE.clear()
        utils.CONFIG_BRIDGE.update(self.config)
        for mod, log in self.loggers:
            mod.log = log

    def translate(self, model=None):
        path = os.path.join(self.tmp.name, "beam.liml")
        with io.open(path, "w", encoding="utf-8") as f:
            f.write((model or beam_model()).xml())
        md = parser.extract_aster_data_from_liml(path)
        self.assertIsNotNone(md, self.logs)
        mail = os.path.join(self.tmp.name, "beam.mail")
        self.assertTrue(mesh.write_mesh(mail, md), self.logs)
        export = writer.generate_AsterSolver_files(
            self.tmp.name, md, mail, "beam")
        return md, export

    def commands(self):
        path = os.path.join(self.tmp.name, utils.run_name(".comm"))
        with io.open(path, encoding="utf-8") as f:
            return f.read()

    def test_full_translation_and_results(self):
        md, export = self.translate()
        self.assertTrue(export, self.logs)
        text = self.commands()
        compile(text, "generated.comm", "exec")
        self.assertIn("MODELISATION='POU_D_T_GD'", text)
        self.assertIn("RELATION='ELAS_POUTRE_GR'", text)
        self.assertIn("DEFORMATION='GROT_GDEP'", text)
        self.assertNotIn("DEFORMATION='PETIT'", text)
        self.assertIn("SIEF_ELNO", text)
        self.assertNotIn("SIPM_ELNO", text)
        self.assertTrue(md["has_rot_dofs"])
        self.assertIn("BLOCAGE=('DEPLACEMENT', 'ROTATION')", text)
        self.assertNotIn("geometrically linear", "\n".join(self.logs))

    def test_general_section_orientation_and_distributed_load(self):
        section = ('<geometric type="UniformGeneralSection" area="0.0006 m2" '
                   'iyy="4.5e-8 m4" izz="2e-8 m4" '
                   'torsionconstant="4e-8 m4" />')
        m = beam_model(section, distributed=True)
        elems = m.elsets[0][2]
        for i, (eid, shape, nids, _) in enumerate(elems):
            elems[i] = (eid, shape, nids, 'orientationvz="1"')
        md, export = self.translate(m)
        self.assertTrue(export, self.logs)
        text = self.commands()
        self.assertIn("'AY', 'AZ'", text)
        self.assertIn("VECT_Y", text)
        self.assertIn("FORCE_POUTRE", text)
        self.assertFalse(md.get("skipped_loads"), md.get("skipped_loads"))

    def test_reject_plastic_material(self):
        _, export = self.translate(beam_model(plastic=True))
        self.assertIsNone(export)
        self.assertIn("requires an elastic material", "\n".join(self.logs))

    def test_reject_linear_or_dynamic_solver(self):
        for analysis in ("lin_static_3d", "dyna_impl_3d"):
            utils.CONFIG_BRIDGE["active_analysis"] = analysis
            _, export = self.translate()
            self.assertIsNone(export)

    def test_reject_small_deformation_and_automatic_imperfection(self):
        a = utils.CONFIG_BRIDGE["analyses"]["nonlin_static_3d"]
        a["deformation"] = "PETIT"
        _, export = self.translate()
        self.assertIsNone(export)
        a["deformation"] = "GROT_GDEP"
        a["imperfection_amplitude"] = 0.001
        _, export = self.translate()
        self.assertIsNone(export)

    def test_classical_beam_regression(self):
        for mod in ("POU_D_E", "POU_D_T"):
            utils.CONFIG_BRIDGE["elem_beam"] = mod
            _, export = self.translate()
            self.assertTrue(export, self.logs)
            text = self.commands()
            self.assertIn("RELATION='ELAS'", text)
            self.assertIn("DEFORMATION='PETIT'", text)
            self.assertIn("SIPM_ELNO", text)
            self.assertNotIn("ELAS_POUTRE_GR", text)

    def test_numerical_benchmark_commands(self):
        import tests_nonlinear_beams as benchmarks
        for build in (benchmarks.small_force,
                      lambda: benchmarks.moment_model(benchmarks.THETA),
                      benchmarks.imperfect_column):
            _, export = self.translate(build())
            self.assertTrue(export, self.logs)
            compile(self.commands(), "benchmark.comm", "exec")

    def test_mixed_beam_and_solid_behaviours(self):
        md, _ = self.translate()
        md["elsets"]["Solid"] = "Steel"
        md["elset_modelisations"]["Solid"] = "3D"
        text = writer._build_comportement_block(md, {"deformation": "GROT_GDEP"})
        self.assertIn("GROUP_MA='Beam', RELATION='ELAS_POUTRE_GR', "
                      "DEFORMATION='GROT_GDEP'", text)
        self.assertIn("GROUP_MA='Solid', RELATION='ELAS', "
                      "DEFORMATION='GREEN_LAGRANGE'", text)
        self.assertEqual(writer._result_groups(md), (["Solid"], [], []))

    def test_gravity_and_line_load_remain_in_load_preparation(self):
        m = beam_model(distributed=True)
        m.add('<gravity><x>0 m/s2</x><y>-9.81 m/s2</y><z>0 m/s2</z></gravity>')
        md, export = self.translate(m)
        self.assertTrue(export, self.logs)
        self.assertIn("FORCE_POUTRE", self.commands())
        self.assertEqual(md["materials"]["Steel"]["rho"], 7850.0)
        self.assertFalse(md.get("skipped_loads"), md.get("skipped_loads"))

    def test_reject_unsupported_element_and_centrifugal_load(self):
        md, _ = self.translate()
        a = utils.CONFIG_BRIDGE["analyses"]["nonlin_static_3d"]
        md["bcs"].append({"aster_keyword": "ROTATION"})
        self.assertFalse(writer._check_gd_beams(md, a))
        md["bcs"].pop()
        eid, _, nodes, group = md["mesh_elements"][0]
        md["mesh_elements"][0] = (eid, "line3", nodes, group)
        self.assertFalse(writer._check_gd_beams(md, a))


if __name__ == "__main__":
    unittest.main()
