"""Compare vertical I-beam orientation with analytical static/modal results.
Usage: python check_beam_orientations.py RUN_ASTER NEW_OUTPUT_DIRECTORY
"""
import json
import math
import subprocess
import sys
from pathlib import Path
from verif_lib import Model, linemesh, Results


def make_model(kind, orientation):
    m = Model('S30' if kind=='static' else 'V30')
    m.elset('Beam','Steel')
    eids,nodes = linemesh(m,'Beam',12,(0,0,0),(0,0,3))
    m.mat('Steel','<geometric type="I" ia="0.1 m" ib="0.05 m" ic="0.01 m" id="0.008 m" />',
          '<mechanical type="Isotropic" youngsmodulus="210 GPa" poissonratio="0.3" />',
          '<density density="7850 kg/m3" />')
    if orientation:
        elems=m.elsets[0][2]
        for i,(eid,shape,nids,_) in enumerate(elems):
            elems[i]=(eid,shape,nids,orientation)
    m.node_selection('Root',[nodes[0]])
    m.node_selection('Tip',[nodes[-1]])
    m.node_selection('Free',nodes[1:])
    m.add('<fixed selection="Root" />')
    # Isolate global-X bending, so a section rotation changes the selected
    # modal family rather than merely reordering two unconstrained modes.
    for tag,unit,dirs in [('displacement','m',((0,1,0),(0,0,1))),
                          ('noderotation','rad',((1,0,0),(0,0,1)))]:
        for x,y,z in dirs:
            m.add('<%s selection="Free"><value>0 %s</value><x>%s</x><y>%s</y><z>%s</z></%s>'%(tag,unit,x,y,z,tag))
    if kind=='static':
        m.add('<force selection="Tip"><x>100 N</x><y>0 N</y><z>0 N</z></force>')
    return m


def main():
    launcher,target=sys.argv[1:3]
    out=Path(target)
    out.mkdir(parents=True,exist_ok=False)
    # Independent rectangular decomposition of the symmetric I section.
    A=2*.05*.01+.08*.008
    Iy=2*.01*.05**3/12+.08*.008**3/12
    Iz=(.05*.1**3-(.05-.008)*.08**3)/12
    report={}
    for name,orientation,I in [('default','',Iy),('vector_y','orientationvy="1"',Iy),
                              ('rotated','orientationvy="1" angle="90 deg"',Iz),
                              ('vector_x','orientationvx="1"',Iz)]:
        for kind in ('static','modal'):
            case=name+'_'+kind
            m=make_model(kind,orientation)
            src=out/(case+'.liml')
            m.write(str(src))
            cfg=out/(case+'.json')
            cfg.write_text(json.dumps({'general':{'elem_beam':'POU_D_E','exe_aster':launcher},
                                      'settings':{'modal_prestress':False,'modal_modes':3},'timeout':120}))
            work=out/(case+'-run')
            assert not work.exists()
            subprocess.run([sys.executable,str(Path(__file__).with_name('_run_one.py')),str(work),str(src),
                            'lin_static_3d' if kind=='static' else 'modal_3d',str(cfg)],check=True)
            run=json.loads((work/'run.json').read_text())
            assert run['diagnostic']=='OK' and not run['errors'],run
            r=Results(run['results'],m)
            actual=r.steps[0][1]['displx'][13] if kind=='static' else r.mode(1)
            expected=100*3**3/(3*210e9*I) if kind=='static' else 1.875104068711961**2/(2*math.pi*3**2)*math.sqrt(210e9*I/(7850*A))
            rel=abs(actual/expected-1)
            report[case]={'actual':actual,'expected':expected,'relative_error':rel}
            print(case,report[case],flush=True)
    (out/'summary.json').write_text(json.dumps(report,indent=2))
    assert all(r['relative_error']<2e-5 for r in report.values()),report
    print('All orientation checks passed.')


if __name__=='__main__':
    main()
