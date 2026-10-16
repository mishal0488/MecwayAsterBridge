"""Offset kinematics, validation and translation regression tests."""
import copy
import os
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT,'Scripts'), os.path.join(ROOT,'TestFiles','Verificatie')]
import bridge_utils as u
import bridge_offsets as offsets
import bridge_from_liml as parser
import bridge_mesh as mesh
import bridge_to_aster as writer
import aster_results_export as results
from check_beam_offsets import model


class OffsetTests(unittest.TestCase):
    def setUp(self):
        self.config = copy.deepcopy(u.CONFIG_BRIDGE)
        u.CONFIG_BRIDGE.update(active_analysis='lin_static_3d', elem_beam='POU_D_E', result_route='liml')
        u.CONFIG_BRIDGE['analyses']['modal_3d']['modal_prestress'] = False
        self.tmp = tempfile.TemporaryDirectory()
        self.logs = []
        self.loggers = [(mod, mod.log) for mod in (u, offsets, parser, mesh, writer)]
        for mod, _ in self.loggers:
            mod.log = self.logs.append

    def tearDown(self):
        self.tmp.cleanup()
        u.CONFIG_BRIDGE.clear()
        u.CONFIG_BRIDGE.update(self.config)
        for mod, logger in self.loggers:
            mod.log = logger

    def parse(self, value=.1, kind='static'):
        path = os.path.join(self.tmp.name,'beam.liml')
        with open(path, 'w') as f:
            f.write(model(kind, value).xml())
        return parser.extract_aster_data_from_liml(path)

    def test_static_and_modal_translation(self):
        for kind, analysis in [('static','lin_static_3d'),('modal','modal_3d')]:
            u.CONFIG_BRIDGE['active_analysis'] = analysis
            md = self.parse(kind=kind)
            mail = os.path.join(self.tmp.name,'beam.mail')
            self.assertTrue(mesh.write_mesh(mail,md), self.logs)
            self.assertEqual(len(md['offset_equations']),120)
            self.assertEqual(len(md['offset_node_aliases']),20)
            self.assertTrue(writer.generate_AsterSolver_files(self.tmp.name,md,mail,'beam'),self.logs)
            with open(os.path.join(self.tmp.name,u.run_name('.comm'))) as f:
                text = f.read()
            compile(text,'generated.comm','exec')
            self.assertIn('LIAISON_DDL',text)
            self.assertIn('offset_nodes=',text)
            self.assertIn("MODELISATION='DIS_TR'",text)

    def test_rigid_motion_for_signed_offsets_and_rotated_axes(self):
        md = self.parse()
        # Rotate the line; all three local offset components must be respected.
        md['mesh_nodes'] = {n:(0.,p[0],0.) for n,p in md['mesh_nodes'].items()}
        md['beam_offsets'][0]['vector'] = (.2,-.3,.4)
        before = dict(md['mesh_nodes'])
        self.assertTrue(offsets.prepare(md,mesh._mecway_v_axis))
        theta, shift = np.array([.03,-.07,.11]), np.array([1.,2.,3.])
        values = {}
        for n,p in md['mesh_nodes'].items():
            vec = np.r_[shift+np.cross(theta,p),theta]
            values.update({(n,d):v for d,v in zip(('DX','DY','DZ','DRX','DRY','DRZ'),vec)})
        for eq in md['offset_equations']:
            self.assertAlmostEqual(sum(values[n,d]*c for n,d,c in eq['terms']),0.,places=12)
        for n,p in before.items():
            self.assertEqual(md['mesh_nodes'][n],p)
        # The same reference junction receives two independent physical nodes.
        aliases = list(md['offset_node_aliases'].values())
        self.assertEqual(aliases.count(2),2)

    def test_global_offset_is_not_rotated_with_beam(self):
        path=os.path.join(self.tmp.name,'global.liml')
        text=model('static',.18).xml().replace('<w>','<z>').replace('</w>','</z>').replace('<u>','<x>').replace('</u>','</x>').replace('<v>','<y>').replace('</v>','</y>')
        with open(path,'w') as f: f.write(text)
        md=parser.extract_aster_data_from_liml(path)
        self.assertEqual(md['beam_offsets'][0]['frame'],'global')
        md['mesh_nodes']={n:(0.,p[0],0.) for n,p in md['mesh_nodes'].items()}
        md['beam_orientations']={e[0]:{'vec':(0.,0.,1.),'angle':0.} for e in md['mesh_elements']}
        before=dict(md['mesh_nodes'])
        self.assertTrue(offsets.prepare(md,mesh._mecway_v_axis))
        for new,old in md['offset_node_aliases'].items():
            np.testing.assert_allclose(np.array(md['mesh_nodes'][new])-before[old],(0,0,.18))
        self.assertEqual(offsets.read_vector(ET.fromstring('<offset><z>180 mm</z></offset>')),(0.,0.,.18))
        self.assertFalse(all(np.isfinite(offsets.read_vector(ET.fromstring('<offset><z>1 m</z><w>2 m</w></offset>')))))

    def test_zero_offsets_preserve_mesh(self):
        md = self.parse(0)
        before = copy.deepcopy(md)
        self.assertTrue(offsets.prepare(md,mesh._mecway_v_axis))
        self.assertEqual(md['mesh_nodes'],before['mesh_nodes'])
        self.assertEqual(md['mesh_elements'],before['mesh_elements'])
        self.assertEqual(md['offset_equations'],[])

    def test_unsupported_modes_and_overlaps_fail(self):
        for analysis in ('nonlin_static_3d','buckling_3d'):
            md = self.parse()
            u.CONFIG_BRIDGE['active_analysis'] = analysis
            # Validation rejects before accessing any analysis-specific fields.
            if analysis in u.CONFIG_BRIDGE['analyses']:
                self.assertFalse(offsets.prepare(md,mesh._mecway_v_axis))
        u.CONFIG_BRIDGE['active_analysis'] = 'modal_3d'
        u.CONFIG_BRIDGE['analyses']['modal_3d']['modal_prestress'] = True
        md=self.parse()
        self.assertFalse(offsets.prepare(md,mesh._mecway_v_axis))
        self.assertIn('disable modal prestress',md['preparation_error'])
        u.CONFIG_BRIDGE['active_analysis'] = 'lin_static_3d'
        md = self.parse()
        md['beam_offsets'] *= 2
        self.assertFalse(offsets.prepare(md,mesh._mecway_v_axis))

    def test_invalid_lengths_do_not_become_zero(self):
        for text in ('bad','1 mm + 2 mm','1 bogus','NaN','1e999 m'):
            vec = offsets.read_vector(ET.fromstring('<offset><w>%s</w></offset>' % text))
            self.assertFalse(np.isfinite(vec[2]))
        self.assertEqual(offsets.read_vector(ET.fromstring('<offset><v>-100 mm</v></offset>')),(0.,-.1,0.))

    def test_suppressed_offset_and_unsupported_load(self):
        path = os.path.join(self.tmp.name,'suppressed.liml')
        with open(path,'w') as f:
            f.write(model('static',.1).xml().replace('<offset ', '<offset suppressed="true" '))
        md = parser.extract_aster_data_from_liml(path)
        self.assertFalse(md.get('beam_offsets'))
        before = copy.deepcopy(md)
        self.assertTrue(offsets.prepare(md,mesh._mecway_v_axis))
        self.assertEqual(md,before)
        md = self.parse()
        md['bcs'].append({'type':'linepressure','selection':'Offset'})
        self.assertFalse(offsets.prepare(md,mesh._mecway_v_axis))

    def test_reactions_sum_multiple_arms_and_moments(self):
        vals = np.zeros((3,6))
        vals[0,0] = -2
        vals[1] = [10,0,0,0,-1,0]
        vals[2] = [20,0,0,0,-4,0]
        out,mask = results._offset_reactions(vals,np.ones_like(vals,dtype=bool),
            ['DX','DY','DZ','DRX','DRY','DRZ'],
            (np.array([1,2]),np.array([0,0]),np.array([[0,0,.1],[0,0,.2]])))
        np.testing.assert_allclose(out[0],[28,0,0,0,0,0],atol=1e-12)
        np.testing.assert_array_equal(vals[0],[-2,0,0,0,0,0])


if __name__ == '__main__':
    unittest.main()
