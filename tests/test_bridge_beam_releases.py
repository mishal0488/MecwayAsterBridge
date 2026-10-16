"""Beam-end release parsing, local constraints and results mapping."""
import copy
import os
import sys
import tempfile
import unittest
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0]=[os.path.join(ROOT,'Scripts'),os.path.join(ROOT,'TestFiles','Verificatie')]
import bridge_utils as u
import bridge_from_liml as parser
import bridge_mesh as mesh
import bridge_to_aster as writer
import bridge_releases as releases
from check_beam_releases import model


class ReleaseTests(unittest.TestCase):
    def setUp(self):
        self.config=copy.deepcopy(u.CONFIG_BRIDGE)
        u.CONFIG_BRIDGE.update(active_analysis='lin_static_3d',elem_beam='POU_D_E',result_route='liml')
        self.tmp=tempfile.TemporaryDirectory()
        self.logs=[]
        self.loggers=[(mod,mod.log) for mod in (u,parser,mesh,writer,releases)]
        for mod,_ in self.loggers: mod.log=self.logs.append

    def tearDown(self):
        self.tmp.cleanup()
        u.CONFIG_BRIDGE.clear()
        u.CONFIG_BRIDGE.update(self.config)
        for mod,logger in self.loggers: mod.log=logger

    def parse(self,axes=('v','w'),text=None):
        path=os.path.join(self.tmp.name,'joint.liml')
        with open(path,'w',encoding='utf-8') as f: f.write(text or model(axes).xml())
        return parser.extract_aster_data_from_liml(path)

    def test_static_and_prestressed_modal_translation(self):
        for analysis in ('lin_static_3d','modal_3d'):
            u.CONFIG_BRIDGE['active_analysis']=analysis
            md=self.parse()
            self.assertFalse(md['skipped_loads'])
            before=dict(md['mesh_nodes'])
            path=os.path.join(self.tmp.name,'joint.mail')
            self.assertTrue(mesh.write_mesh(path,md),self.logs)
            self.assertEqual(md['offset_node_aliases'],{22:11})
            self.assertEqual(md['mesh_nodes'][22],before[11])
            self.assertEqual(len(md['release_equations']),4)
            self.assertTrue(writer.generate_AsterSolver_files(self.tmp.name,md,path,'joint'),self.logs)
            with open(os.path.join(self.tmp.name,u.run_name('.comm'))) as f: text=f.read()
            compile(text,'generated.comm','exec')
            self.assertIn('LIAISON_DDL',text)
            self.assertIn('offset_nodes={22: 11}',text)
            # No reference-node springs or penalty stiffness are introduced.
            self.assertNotIn("MODELISATION='DIS_TR'",text)

    def test_rotated_single_axis_releases_correct_local_rotation(self):
        md=self.parse(('v',))
        md['beam_orientations'][11]={'vec':(0.,1.,1.),'angle':0.}
        self.assertTrue(releases.prepare(md,mesh._mecway_v_axis))
        for direction,expected_free in [((0.,1.,1.),True),((1.,0.,0.),False),((0.,-1.,1.),False)]:
            vals=dict(((22,d),v) for d,v in zip(('DRX','DRY','DRZ'),direction))
            residuals=[sum(vals.get((n,d),0.)*c for n,d,c in eq['terms']) for eq in md['release_equations']]
            self.assertEqual(all(abs(v)<1e-12 for v in residuals),expected_free)
        self.assertEqual(len(md['release_equations']),5)

    def test_all_rotations_released_preserves_translations(self):
        md=self.parse(('u','v','w'))
        self.assertTrue(releases.prepare(md,mesh._mecway_v_axis))
        self.assertEqual(len(md['release_equations']),3)
        self.assertTrue(all(d in ('DX','DY','DZ') for eq in md['release_equations'] for n,d,c in eq['terms']))

    def test_suppression_and_no_release_are_noops(self):
        text=model().xml().replace('<flexiblejoint ', '<flexiblejoint suppressed="true" ')
        md=self.parse(text=text)
        before=copy.deepcopy(md)
        self.assertTrue(releases.prepare(md,mesh._mecway_v_axis))
        self.assertEqual(md,before)

    def test_invalid_faces_and_unsupported_cases_stop_without_mesh_changes(self):
        cases=[('face',None),('isolated',None),('offset',None),('analysis','nonlin_static_3d')]
        for case,analysis in cases:
            u.CONFIG_BRIDGE['active_analysis']='lin_static_3d'
            md=self.parse()
            if case=='face': md['face_selection_faces']['Release']=[(11,1)]
            if case=='isolated': md['face_selection_faces']['Release']=[(1,2)]
            if case=='offset': md['beam_offsets']=[{'selection':'Beam','vector':(0,.1,0)}]
            if analysis: u.CONFIG_BRIDGE['active_analysis']=analysis
            before=copy.deepcopy(md['mesh_nodes']),copy.deepcopy(md['mesh_elements'])
            self.assertFalse(releases.prepare(md,mesh._mecway_v_axis),case)
            self.assertEqual((md['mesh_nodes'],md['mesh_elements']),before)

    def test_unsupported_joint_parameters_not_silently_ignored(self):
        text=model().xml().replace('<rotv />','<rotv>123</rotv>')
        self.assertIsNone(self.parse(text=text))


if __name__=='__main__': unittest.main()
