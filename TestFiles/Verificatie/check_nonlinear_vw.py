"""Native end-to-end checks for exact zero-twist nonlinear V/W joints.

Usage: python check_nonlinear_vw.py RUN_ASTER NEW_OUTPUT_DIRECTORY [CASE ...]
Cases: compound oblique refined torsion both_ends both_ends_compound small linear_small
"""
import json, math, subprocess, sys
from pathlib import Path
import numpy as np
from verif_lib import Model, Results, linemesh
import check_nonlinear_releases as sphere


def model(kind):
    if kind in ('both_ends','both_ends_compound'):
        m=Model('L30');m.elset('Beam','Steel')
        es,ns=linemesh(m,'Beam',18,(0,0,0),(3,0,0))
        m.mat('Steel','<geometric type="UniformGeneralSection" area="0.0001 m2" iyy="%s m4" izz="%s m4" torsionconstant="%s m4" />'%(sphere.I,sphere.I,sphere.J),
            '<mechanical type="Isotropic" youngsmodulus="210 GPa" poissonratio="0.3" />')
        m.node_selection('Roots',[ns[0],ns[-1]])
        m.node_selection('Tip',[ns[9]])
        m.face_selection('Release',[(es[6],2),(es[11],3)])
        m.add('<fixed selection="Roots" />')
        m.add('<flexiblejoint selection="Release"><rotv/><rotw/></flexiblejoint>')
        torque=sphere.G*sphere.J*math.pi/3
        m.add('<moment selection="Tip"><x time="0 s">0 N.m</x><x time="1 s">%.15g N.m</x><y>0 N.m</y><z>0 N.m</z></moment>'%torque)
        if kind=='both_ends_compound':
            m.add('<force selection="Tip"><x>0 N</x><y time="0 s">0 N</y><y time="1 s">700 N</y><z time="0 s">0 N</z><z time="1 s">90 N</z></force>')
        return m,ns[9]
    sphere.N=24 if kind in ('small','linear_small') else 6
    m,Q,roots,load=sphere.model('oblique' if kind=='oblique' else 'torsion' if kind=='torsion' else 'bending')
    m.items=[s.replace('<rotu/>','') for s in m.items]
    if kind!='torsion':
        factor=1e-4 if kind in ('small','linear_small') else 1.
        if factor!=1.:
            m.items=[s.replace('700 N','0.07 N') for s in m.items]
        moment=Q@np.array([50.,25.,70.])*factor
        m.add('<moment selection="Tip">'+''.join('<%s time="0 s">0 N.m</%s><%s time="1 s">%.15g N.m</%s>'%(c,c,c,v,c) for c,v in zip('xyz',moment))+'</moment>')
    return m,sphere.N+1


def main():
    launcher,target=sys.argv[1:3]
    kinds=sys.argv[3:] or ['compound','oblique','refined','torsion','both_ends','both_ends_compound','small','linear_small']
    out=Path(target).resolve();out.mkdir(parents=True,exist_ok=False)
    report={}
    for kind in kinds:
        m,tip=model(kind);src=out/(kind+'.liml');m.write(str(src))
        steps=80 if kind=='refined' else 40
        linear=kind=='linear_small'
        cfg=out/(kind+'.json')
        cfg.write_text(json.dumps({'general':{'elem_beam':'POU_D_T' if linear else 'POU_D_T_GD','exe_aster':launcher},
            'settings':{'deformation':'GROT_GDEP','substeps_per_interval':steps,'newton_max_iter':50},'timeout':900}))
        work=out/(kind+'-run');assert not work.exists()
        subprocess.run([sys.executable,str(Path(__file__).with_name('_run_one.py')),str(work),str(src),'lin_static_3d' if linear else 'nonlin_static_3d',str(cfg)],check=True)
        run=json.loads((work/'run.json').read_text())
        assert run['results'] and not run['errors'],run
        r=Results(run['results'],m)
        t,d=r.steps[-1]
        values=[d[field][tip] for field in ('displx','disply','displz','rotx','roty','rotz')]
        if not linear:
            assert abs(t-1)<1e-9 and len(r.steps)==steps,(t,len(r.steps))
        if kind in ('torsion','both_ends'):
            np.testing.assert_allclose(values[3],math.pi/4,rtol=1e-7)
        report[kind]={'tip':values,'steps':len(r.steps),'seconds':run['seconds']}
        (out/'summary.json').write_text(json.dumps(report,indent=2))
        print(kind,report[kind],flush=True)
    if 'small' in report and 'linear_small' in report:
        np.testing.assert_allclose(report['small']['tip'],report['linear_small']['tip'],rtol=.001,atol=2e-8)
    if 'compound' in report and 'oblique' in report:
        Q=sphere.model('oblique')[1]
        for i in (0,3):
            np.testing.assert_allclose(Q.T@np.array(report['oblique']['tip'][i:i+3]),report['compound']['tip'][i:i+3],atol=1e-8)
    print('Requested nonlinear V/W end-to-end checks passed.')

if __name__=='__main__':main()
