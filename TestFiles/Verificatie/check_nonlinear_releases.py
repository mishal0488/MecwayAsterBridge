"""Exact nonlinear spherical beam-joint checks (no offsets or penalty springs).
Usage: python check_nonlinear_releases.py RUN_ASTER NEW_OUTPUT_DIRECTORY
Requires NumPy. Selective releases are intentionally not enabled.
"""
import json, math, subprocess, sys
from pathlib import Path
import numpy as np
from verif_lib import Model, Results

E=210.e9
G=E/2.6
I=1.e-8/12
J=1.4e-9
N=40


def model(kind):
    m=Model('L30')
    m.elset('Beam','Steel')
    m.mat('Steel','<geometric type="UniformGeneralSection" area="0.0001 m2" iyy="%s m4" izz="%s m4" torsionconstant="%s m4" />'%(I,I,J),
          '<mechanical type="Isotropic" youngsmodulus="210 GPa" poissonratio="0.3" />')
    # Two independent, coincident cantilevers meet only at their tips.
    # A spherical release on the second tip must transmit force but no torque.
    Q=np.eye(3)
    if kind=='oblique':
        a=.7;b=.6
        Q=np.array([[math.cos(a),-math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,1]]) @ np.array([[math.cos(b),0,math.sin(b)],[0,1,0],[-math.sin(b),0,math.cos(b)]])
    roots=[]
    for beam in range(2):
        ids=[]
        for i in range(N+1):
            if beam and i==N: ids.append(N+1);continue
            m.nodes.append(tuple(Q @ np.array([i/N,0.,0.])))
            ids.append(len(m.nodes))
        roots.append(ids[0])
        for i in range(N):m.elem('Beam','line2',ids[i:i+2])
    m.node_selection('Roots',roots)
    m.node_selection('Tip',[N+1])
    m.face_selection('Release',[(2*N,3)])
    m.add('<fixed selection="Roots" />')
    m.add('<flexiblejoint selection="Release"><rotu/><rotv/><rotw/></flexiblejoint>')
    load=Q @ np.array([G*J*math.pi/2,0,0] if kind=='torsion' else [0,700.,0])
    tag,unit=('moment','N.m') if kind=='torsion' else ('force','N')
    m.add('<%s selection="Tip">'%tag+''.join('<%s time="0 s">0 %s</%s><%s time="1 s">%.14g %s</%s>'%(c,unit,c,c,v,unit,c) for c,v in zip('xyz',load))+'</%s>'%tag)
    return m,Q,roots,load


def main():
    launcher,target=sys.argv[1:3]
    out=Path(target).resolve();out.mkdir(parents=True,exist_ok=False)
    # Independent planar inextensible elastica: theta''=-F*cos(theta)/EI,
    # theta(0)=0, theta'(1)=0, x'=cos(theta), y'=sin(theta).
    def shoot(curvature):
        y=np.array([0.,curvature,0.,0.]);h=1/1200
        def rhs(y):return np.array([y[1],-350/(E*I)*math.cos(y[0]),math.cos(y[0]),math.sin(y[0])])
        for _ in range(1200):
            a=rhs(y);b=rhs(y+h*a/2);c=rhs(y+h*b/2);d=rhs(y+h*c)
            y+=h*(a+2*b+2*c+d)/6
        return y
    low,high=0.,350/(E*I)
    assert shoot(low)[1]<0<shoot(high)[1]
    for _ in range(45):
        mid=(low+high)/2
        if shoot(mid)[1]>0:high=mid
        else:low=mid
    solved=shoot((low+high)/2)
    reference=np.array([solved[2]-1,solved[3],0])
    cfg=out/'settings.json'
    cfg.write_text(json.dumps({'general':{'elem_beam':'POU_D_T_GD','exe_aster':launcher},'settings':{'deformation':'GROT_GDEP','substeps_per_interval':40,'newton_max_iter':50,'convergence_resi_rela':1e-7,'convergence_criterion':'RELA_OF_MAXI','convergence_resi_maxi':1e-5},'timeout':180}))
    report={}
    for kind in ('bending','oblique','torsion'):
        m,Q,roots,load=model(kind)
        src=out/(kind+'.liml');m.write(str(src));work=out/(kind+'-run')
        assert not work.exists()
        subprocess.run([sys.executable,str(Path(__file__).with_name('_run_one.py')),str(work),str(src),'nonlin_static_3d',str(cfg)],check=True)
        run=json.loads((work/'run.json').read_text())
        assert run['diagnostic'] in ('OK','<A>_ALARM') and not run['errors'] and run['results'],run
        r=Results(run['results'],m);t,d=r.steps[-1];assert abs(t-1)<1e-10
        tip=np.array([d['displ'+c][N+1] for c in 'xyz'])
        if kind=='torsion':
            np.testing.assert_allclose(d['rotx'][N+1],math.pi/2,rtol=1e-5)
            np.testing.assert_allclose(d['reactionmomx'][roots[0]],-load[0],rtol=1e-6)
            np.testing.assert_allclose(d['reactionmomx'][roots[1]],0,atol=1e-6)
        else:
            # Elastica is inextensible and shear-rigid; these corrections are
            # tiny for this slender beam. Allow 0.15% for mesh + shear effects.
            np.testing.assert_allclose(Q.T @ tip,reference,rtol=.0015,atol=2e-5)
            total=np.array([sum(d['reactionforce'+c][n] for n in roots) for c in 'xyz'])
            np.testing.assert_allclose(total,-load,atol=1e-5)
        # Inspect the generated end node directly. POU_D_T_GD SIEF_ELNO
        # extrapolates integration-point section moments; those are not a
        # reliable nodal hinge-equilibrium check on a finite element mesh.
        comm=Path(run['export']).with_suffix('.comm')
        raw=work/'joint_raw.json'
        extra="import json\n_raw={}\n"
        for field in ('DEPL','REAC_NODA'):
            extra+="_sf=res.getField(%r,max(res.getIndexes())).toSimpleFieldOnNodes()\n"%field
            extra+="_vals,_mask=_sf.getValues()\n"
            extra+="_raw[%r]=dict(zip(_sf.getComponents(),[float(v) for v in _vals[%d]]))\n"%(field,2*N+1)
        extra+="with open(%r,'w') as _f: json.dump(_raw,_f)\n"%str(raw)
        comm.write_text(comm.read_text().replace('FIN()',extra+'FIN()'))
        with (work/'raw_console.txt').open('w') as f:
            subprocess.run([launcher,run['export']],stdout=f,stderr=subprocess.STDOUT,timeout=180)
        assert raw.exists(), 'Raw joint export did not complete'
        joint=json.loads(raw.read_text())
        for c in ('DRX','DRY','DRZ'):
            np.testing.assert_allclose(joint['REAC_NODA'][c],0,atol=1e-5)
        for c in 'xyz':
            np.testing.assert_allclose(joint['DEPL']['D'+c.upper()],d['displ'+c][N+1],atol=1e-9)
        if kind=='torsion':
            np.testing.assert_allclose(joint['DEPL']['DRX'],0,atol=1e-7)
        report[kind]={'tip':tip.tolist(),'reference_elastica':(Q@reference).tolist() if kind!='torsion' else None,'rotation_x':d['rotx'][N+1]}
        print(kind,report[kind],flush=True)
    np.testing.assert_allclose(report['bending']['tip'],model('oblique')[1].T @ np.array(report['oblique']['tip']),atol=1e-7)
    (out/'summary.json').write_text(json.dumps(report,indent=2))
    print('All exact nonlinear spherical-joint checks passed.')

if __name__=='__main__':main()
