"""Numerical offset benchmarks. Usage: python check_beam_offsets.py LAUNCHER OUTPUT_DIR.
OUTPUT_DIR must not exist. Requires numpy and a working Code_Aster launcher.
"""
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np

from verif_lib import Model, linemesh, Results, _values


def model(kind, offset):
    m = Model('S30' if kind == 'static' else 'V30')
    m.elset('Beam', 'Steel')
    eids, nodes = linemesh(m, 'Beam', 10, (0, 0, 0), (1, 0, 0))
    m.mat('Steel', '<geometric type="RectBar" recta="0.02 m" rectb="0.03 m" />',
          '<mechanical type="Isotropic" youngsmodulus="210 GPa" poissonratio="0.3" />',
          '<density density="7850 kg/m3" />')
    m.element_selection('Offset', eids)
    m.node_selection('Root', [nodes[0]])
    m.node_selection('Tip', [nodes[-1]])
    m.node_selection('Free', nodes[1:])
    m.add('<offset selection="Offset"><u>0 m</u><v>0 m</v><w>%s m</w></offset>' % offset)
    m.add('<fixed selection="Root" />')
    if kind == 'static':
        m.add('<force selection="Tip"><x>100 N</x><y>0 N</y><z>0 N</z></force>')
    if kind == 'restrained':
        # Only reference-node DRY is free: eccentric axial mass and stiffness
        # both contribute. This is deliberately not a translated cantilever.
        for tag, unit, directions in [('displacement', 'm', ((1,0,0),(0,1,0),(0,0,1))),
                                      ('noderotation', 'rad', ((1,0,0),(0,0,1)))]:
            for x,y,z in directions:
                m.add('<%s selection="Free"><value>0 %s</value><x>%s</x><y>%s</y><z>%s</z></%s>' % (tag,unit,x,y,z,tag))
    return m


def reference(offset):
    # Independent Euler-Bernoulli consistent stiffness/mass assembly, reduced
    # analytically by ux_physical = offset * theta_y and uz_physical = 0.
    E, A, I, rho, length = 210e9, .02*.03, .02*.03**3/12, 7850., .1
    K, M = np.zeros((11,11)), np.zeros((11,11))
    ke = E*A*offset**2/length*np.array([[1.,-1.],[-1.,1.]]) + E*I/length*np.array([[4.,2.],[2.,4.]])
    me = rho*A*offset**2*length/6*np.array([[2.,1.],[1.,2.]]) + rho*A*length**3/420*np.array([[4.,-3.],[-3.,4.]])
    for i in range(10):
        K[i:i+2,i:i+2] += ke
        M[i:i+2,i:i+2] += me
    invL = np.linalg.inv(np.linalg.cholesky(M[1:,1:]))
    return np.sqrt(np.linalg.eigvalsh(invL @ K[1:,1:] @ invL.T))[:6]/(2*np.pi)


def main():
    launcher, target = sys.argv[1:3]
    out = Path(target)
    out.mkdir(parents=True, exist_ok=False)
    report = {}
    for name, kind, offset in [('static','static',.1), ('modal','modal',.1),
                               ('modal_zero','modal',0), ('restrained','restrained',.1),
                               ('static_t','static',.1), ('modal_t','modal',.1),
                               ('modal_t_zero','modal',0)]:
        m = model(kind, offset)
        src = out/(name+'.liml')
        m.write(str(src))
        cfg = out/(name+'.json')
        beam = 'POU_D_T' if name in ('static_t','modal_t','modal_t_zero') else 'POU_D_E'
        cfg.write_text(json.dumps({'general':{'elem_beam':beam,'exe_aster':launcher},
                                  'settings':{'modal_prestress':False,'modal_modes':6},'timeout':120}))
        work = out/(name+'-run')
        assert not work.exists()
        subprocess.run([sys.executable,str(Path(__file__).with_name('_run_one.py')),str(work),str(src),
                        'lin_static_3d' if kind=='static' else 'modal_3d',str(cfg)],check=True)
        run = json.loads((work/'run.json').read_text())
        assert run['diagnostic']=='OK' and not run['errors'], run
        r = Results(run['results'], m)
        if kind == 'static':
            data = r.steps[0][1]
            expected = 100/(210e9*.0006)+100*.1**2/(210e9*4.5e-8)
            np.testing.assert_allclose(data['displx'][11], expected, rtol=1e-8)
            np.testing.assert_allclose(data['reactionforcex'][1], -100., rtol=1e-8)
            for key in ('reactionforcemag','reactionmommag'):
                np.testing.assert_allclose(data[key][11], 0., atol=1e-7)
            np.testing.assert_allclose(data['reactionmommag'][1], 0., atol=1e-7)
            fields = {v.get('name'): _values(v) for v in
                      ET.parse(run['results']).getroot().iter('elementnodevariable')}
            for key, value in [('tensileforce',100.),('bendingmomv',-10.),
                               ('sheary',0.),('shearz',0.),('torsion',0.),('bendingmomw',0.)]:
                assert len(fields[key]) == 20, (key, len(fields[key]))
                np.testing.assert_allclose(fields[key], value, rtol=1e-8, atol=1e-7)
            report[name] = {k: {str(n):v for n,v in vals.items() if n in (1,11)} for k,vals in data.items() if 'react' in k or k=='displx'}
        else:
            report[name] = [r.modes[i] for i in range(1,7)]
        print(name, report[name], flush=True)
    np.testing.assert_allclose(report['modal'], report['modal_zero'], rtol=1e-7)
    np.testing.assert_allclose(report['modal_t'], report['modal_t_zero'], rtol=1e-7)
    expected = reference(.1)
    report['restrained_reference'] = expected.tolist()
    print('Independent modal reference:', expected, flush=True)
    np.testing.assert_allclose(report['restrained'], expected, rtol=1e-5)
    (out/'summary.json').write_text(json.dumps(report, indent=2))
    print('All numerical offset checks passed.')


if __name__ == '__main__':
    main()
