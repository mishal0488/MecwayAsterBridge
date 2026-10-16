"""Independent beam-end release checks; run with RUN_ASTER and a NEW output directory."""
import json
import math
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from verif_lib import Model, linemesh, Results, _values


def model(axes=('v','w'), modal=False, load='y'):
    m=Model('V30' if modal else 'S30')
    m.elset('Beam','Steel')
    left,ln=linemesh(m,'Beam',10,(0,0,0),(1,0,0))
    right,rn=linemesh(m,'Beam',10,(1,0,0),(2,0,0))
    m.mat('Steel','<geometric type="UniformGeneralSection" area="0.0006 m2" iyy="4.5e-8 m4" izz="2e-8 m4" torsionconstant="4e-8 m4" />',
          '<mechanical type="Isotropic" youngsmodulus="210 GPa" poissonratio="0.3" />',
          '<density density="7850 kg/m3" />')
    m.node_selection('Supports',[ln[0],rn[-1]])
    m.node_selection('Joint',[ln[-1]])
    m.face_selection('Release',[(right[0],2)])
    m.add('<fixed selection="Supports" />')
    m.add('<flexiblejoint selection="Release">'+''.join('<rot%s />'%a for a in axes)+'</flexiblejoint>')
    if not modal:
        if load=='torsion':
            m.add('<moment selection="Joint"><x>100 N.m</x><y>0 N.m</y><z>0 N.m</z></moment>')
        else:
            m.add('<force selection="Joint"><x>0 N</x><y>%s N</y><z>%s N</z></force>' % (100 if load=='y' else 0,100 if load=='z' else 0))
    return m


def main():
    launcher,target=sys.argv[1:3]
    out=Path(target)
    out.mkdir(parents=True,exist_ok=False)
    report={}
    cases=[('hinge',('v','w'),False,'y','POU_D_E'),
           ('single_v',('v',),False,'z','POU_D_E'),
           ('single_w',('w',),False,'y','POU_D_E'),
           ('torsion',('u',),False,'torsion','POU_D_E'),
           ('ball',('u','v','w'),False,'y','POU_D_E'),
           ('modal',('v','w'),True,'y','POU_D_E'),
           ('hinge_t',('v','w'),False,'y','POU_D_T')]
    for name,axes,modal,load,beam in cases:
        m=model(axes,modal,load)
        src=out/(name+'.liml')
        m.write(str(src))
        cfg=out/(name+'.json')
        cfg.write_text(json.dumps({'general':{'elem_beam':beam,'exe_aster':launcher},
                                  'settings':{'modal_prestress':False,'modal_modes':4},'timeout':120}))
        work=out/(name+'-run')
        assert not work.exists()
        subprocess.run([sys.executable,str(Path(__file__).with_name('_run_one.py')),str(work),str(src),
                        'modal_3d' if modal else 'lin_static_3d',str(cfg)],check=True)
        run=json.loads((work/'run.json').read_text())
        assert run['diagnostic']=='OK' and not run['errors'],run
        r=Results(run['results'],m)
        if modal:
            actual=r.mode(1)
            expected=1.875104068711961**2/(2*math.pi)*math.sqrt(210e9*2e-8/(7850*.0006))
            np.testing.assert_allclose(actual,expected,rtol=2e-5)
        else:
            d=r.steps[0][1]
            expected=100/(210e9/2.6*4e-8) if load=='torsion' else 100/(6*210e9*(2e-8 if load=='y' else 4.5e-8))
            if beam=='POU_D_T': expected+=100/(2*(210e9/2.6)*.0006)
            actual=d['rotx' if load=='torsion' else 'displ'+load][11]
            np.testing.assert_allclose(actual,expected,rtol=1e-8)
            field= 'torsion' if load=='torsion' else ('bendingmomw' if load=='y' else 'bendingmomv')
            vals=_values(ET.parse(run['results']).getroot().find(".//elementnodevariable[@name='%s']"%field))
            # Right beam's first end is the released one (element 11).
            np.testing.assert_allclose(vals[20],0.,atol=1e-7)
            component='reactionmomx' if load=='torsion' else 'reactionforce'+load
            np.testing.assert_allclose(sum(d[component].values()),-100.,atol=1e-7)
        report[name]={'actual':actual,'expected':expected}
        print(name,report[name],flush=True)
    (out/'summary.json').write_text(json.dumps(report,indent=2))
    print('All beam-end release checks passed.')


if __name__=='__main__': main()
