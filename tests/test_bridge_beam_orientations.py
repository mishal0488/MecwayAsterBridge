"""Every beam carries Mecway local axes, independent of component contents."""
import copy
import math
import os
import re
import sys
import tempfile
import unittest
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0]=[os.path.join(ROOT,'Scripts'),os.path.join(ROOT,'TestFiles','Verificatie')]
import bridge_utils as u
import bridge_from_liml as parser
import bridge_mesh as mesh
import bridge_to_aster as writer
from verif_lib import Model


class OrientationTests(unittest.TestCase):
    def setUp(self):
        self.config=copy.deepcopy(u.CONFIG_BRIDGE)
        u.CONFIG_BRIDGE.update(active_analysis='lin_static_3d',elem_beam='POU_D_E')
        self.tmp=tempfile.TemporaryDirectory()
        self.logs=[]
        self.loggers=[(mod,mod.log) for mod in (u,parser,mesh,writer)]
        for mod,_ in self.loggers: mod.log=self.logs.append

    def tearDown(self):
        self.tmp.cleanup()
        u.CONFIG_BRIDGE.clear()
        u.CONFIG_BRIDGE.update(self.config)
        for mod,logger in self.loggers: mod.log=logger

    def translate(self, beams):
        m=Model('S30')
        m.mat('Steel','<geometric type="RectBar" recta="0.02 m" rectb="0.03 m" />',
              '<mechanical type="Isotropic" youngsmodulus="210 GPa" poissonratio="0.3" />')
        for group,p,q,attrs in beams:
            m.elset(group,'Steel')
            m.elem(group,'line2',[m.node(*p),m.node(*q)],attrs)
        m.node_selection('Root',[1])
        m.add('<fixed selection="Root" />')
        m.node_selection('Tip',[2])
        m.add('<force selection="Tip"><x>0 N</x><y>1 N</y><z>0 N</z></force>')
        src=os.path.join(self.tmp.name,'orient.liml')
        with open(src,'w') as f: f.write(m.xml())
        md=parser.extract_aster_data_from_liml(src)
        mail=os.path.join(self.tmp.name,'orient.mail')
        self.assertTrue(mesh.write_mesh(mail,md),self.logs)
        self.assertTrue(writer.generate_AsterSolver_files(self.tmp.name,md,mail,'orient'),self.logs)
        with open(os.path.join(self.tmp.name,u.run_name('.comm'))) as f: text=f.read()
        compile(text,'generated.comm','exec')
        return md,text

    def test_defaults_all_axes_and_components(self):
        cases=[('X',(1,0,0),(0,1,0)),('Y',(0,1,0),(-1,0,0)),
               ('Z',(0,0,1),(0,1,0)),('MinusZ',(0,0,-1),(0,1,0)),
               ('Slope',(1,1,1),(-1/math.sqrt(2),1/math.sqrt(2),0))]
        md,text=self.translate([(name,(0,0,0),p,'') for name,p,v in cases])
        self.assertEqual(len(md['beam_vect_y']),len(cases))
        expected={name:v for name,p,v in cases}
        for group,v in zip(sorted(expected),md['beam_vect_y'].values()):
            for a,b in zip(v,expected[group]): self.assertAlmostEqual(a,b,places=12)
        self.assertEqual(len(re.findall(r"_F\(GROUP_MA='BOR\d+', CARA='VECT_Y'",text)),len(cases))

    def test_mixed_default_and_explicit_components(self):
        md,text=self.translate([('Default',(0,0,0),(1,0,0),''),
                                ('Explicit',(0,0,0),(0,1,0),'orientationvz="1"'),
                                ('Rotated',(0,0,0),(0,0,1),'orientationvy="1" angle="90 deg"')])
        for actual,expected in zip(md['beam_vect_y'].values(),[(0,1,0),(0,0,1),(-1,0,0)]):
            for a,b in zip(actual,expected): self.assertAlmostEqual(a,b,places=12)
        self.assertEqual(len(md['beam_vect_y']),3)

    def test_nonbeam_line_elements_get_no_orientation(self):
        md,text=self.translate([('Beam',(0,0,0),(1,0,0),''),
                                ('Truss',(0,1,0),(1,1,0),'truss="true" orientationvz="1"')])
        self.assertEqual(len(md['beam_vect_y']),1)
        self.assertEqual(len(re.findall(r"_F\(GROUP_MA='BOR\d+', CARA='VECT_Y'",text)),1)

    def test_geometric_nonlinear_beams_also_receive_defaults(self):
        u.CONFIG_BRIDGE.update(active_analysis='nonlin_static_3d',elem_beam='POU_D_T_GD')
        md,text=self.translate([('Beam',(0,0,0),(0,0,1),'')])
        self.assertEqual(list(md['beam_vect_y'].values()),[(0.,1.,0.)])
        self.assertIn("MODELISATION='POU_D_T_GD'",text)


if __name__=='__main__': unittest.main()
