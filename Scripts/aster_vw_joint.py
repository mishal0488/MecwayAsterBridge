# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Exact zero-twist V/W joints for elastic POU_D_T_GD static analysis.

Runs inside Code_Aster's Python, not Mecway's IronPython. Each accepted load
increment satisfies the quaternion constraint and its spatial moment gradient.
Native Lagrange multipliers impose the trial constraints; no penalty stiffness
or connection length is introduced. Intermediate outer iterations are discarded.
"""
import numpy as np
from code_aster import CA
from code_aster.Commands import (AFFE_CHAR_MECA_F, DEFI_FONCTION, STAT_NON_LINE,
                                 CALC_CHAMP, DETRUIRE, DEFI_GROUP)
from code_aster.Cata.Syntax import _F
from aster_joint_kinematics import from_vector, constraint


def solve(joints, times, **kwargs):
    times = np.asarray(times, dtype=float)
    if len(times) < 2 or not np.all(np.diff(times) > 0):
        raise ValueError('V/W joints require increasing load times.')
    mesh = kwargs['MODELE'].getMesh()
    nodes = sorted(set(n for j in joints for n in (j['reference'], j['end'])))
    groups = ['__VW%06d'%i for i in range(len(nodes))]
    if set(groups).intersection(mesh.getGroupsOfNodes()):
        raise ValueError('Reserved __VW node groups already exist.')
    DEFI_GROUP(reuse=mesh,MAILLAGE=mesh,CREA_GROUP_NO=tuple(
        _F(NOM=g,NOEUD='N%d'%n) for g,n in zip(groups,nodes)))
    rows = [mesh.getNodes(g)[0] for g in groups]
    indices = {n: i for i, n in enumerate(nodes)}
    joint_data = [(indices[j['reference']], indices[j['end']], np.asarray(j['axis'], dtype=float)) for j in joints]
    previous = None
    retained = []
    x = np.zeros(3 * len(nodes))
    result = CA.NonLinearResult()
    result.allocate(len(times)-1)
    excit = kwargs['EXCIT']
    base_excit = tuple(excit) if isinstance(excit,(list,tuple)) else (excit,)

    def kinematics(q, a, b, axis):
        return constraint(from_vector(q[3*a:3*a+3]), from_vector(q[3*b:3*b+3]), axis)

    def read(r, field, idx):
        sf = r.getField(field, idx).toSimpleFieldOnNodes()
        vals, mask = sf.getValues()
        cols = [list(sf.getComponents()).index(c) for c in ('DRX', 'DRY', 'DRZ')]
        selection = np.ix_(rows, cols)
        if not np.all(mask[selection]):
            raise RuntimeError('V/W joint rotational fields are missing.')
        return np.asarray(vals[selection]).reshape(-1)

    for stage, (start, end) in enumerate(zip(times[:-1], times[1:]), 1):
        old = x.copy()
        inverse = np.eye(len(x))
        last_x = last_f = None
        for iteration in range(40):
            functions, equations = [], []
            for a, b, axis in joint_data:
                c, g = kinematics(x, a, b, axis)
                subset = np.r_[np.arange(3*a,3*a+3), np.arange(3*b,3*b+3)]
                fn = DEFI_FONCTION(NOM_PARA='INST',
                    VALE=(float(start),float(g @ old[subset]),float(end),float(g @ x[subset]-c)),
                    PROL_DROITE='CONSTANT',PROL_GAUCHE='CONSTANT')
                functions.append(fn)
                equations.append(_F(NOEUD=tuple(['N%d'%nodes[a]]*3+['N%d'%nodes[b]]*3),
                    DDL=('DRX','DRY','DRZ','DRX','DRY','DRZ'),
                    COEF_MULT=tuple(float(v) for v in g),COEF_IMPO=fn))
            bc = AFFE_CHAR_MECA_F(MODELE=kwargs['MODELE'],LIAISON_DDL=tuple(equations))
            args = dict(kwargs)
            args['EXCIT'] = base_excit + (_F(CHARGE=bc),)
            args['INCREMENT'] = _F(**dict(kwargs['INCREMENT'],INST_INIT=float(start),INST_FIN=float(end)))
            if previous is not None:
                args['ETAT_INIT'] = _F(EVOL_NOLI=previous,INST=float(start))
            # Never recover partial results after an inner or outer failure.
            r = STAT_NON_LINE(**args)
            idx = max(r.getIndexes())
            if abs(r.getAccessParameters()['INST'][-1]-end) > 1e-10*max(1.,abs(end)):
                raise RuntimeError('V/W trial did not reach its requested load time.')
            actual = read(r, 'DEPL', idx)
            f = actual-x
            closure = max(abs(kinematics(actual,a,b,u)[0]) for a,b,u in joint_data)
            print('BRIDGE_VW_ITER',stage,iteration,float(max(abs(f))),float(closure),flush=True)
            if max(abs(f)) < 1e-8 and closure < 1e-9:
                r = CALC_CHAMP(reuse=r,RESULTAT=r,FORCE='REAC_NODA')
                moments = read(r,'REAC_NODA',idx)
                for a,b,u in joint_data:
                    _,g = kinematics(actual,a,b,u)
                    moment = moments[3*b:3*b+3]
                    error = np.linalg.norm(np.cross(moment,g[3:]))/max(1.,np.linalg.norm(g[3:]))
                    if error > 1e-5 + 1e-7*np.linalg.norm(moment):
                        raise RuntimeError('V/W joint moment equilibrium failed: %g' % error)
                previous = r
                retained.append((r,bc)+tuple(functions))
                x = actual
                # Keep only states satisfying the custom constraint, and retain
                # their own loads for subsequent native reaction recovery.
                result.setModel(kwargs['MODELE'],stage)
                result.setMaterialField(kwargs['CHAM_MATER'],stage)
                result.setElementaryCharacteristics(kwargs['CARA_ELEM'],stage)
                result.setListOfLoads(r.getListOfLoads(idx),stage)
                result.setTime(float(end),stage)
                for field in ('DEPL','SIEF_ELGA','VARI_ELGA','COMPORTEMENT','REAC_NODA'):
                    result.setField(r.getField(field,idx),field,stage)
                print('BRIDGE_VW_CONVERGED',stage,float(end),flush=True)
                break
            DETRUIRE(NOM=(r,bc)+tuple(functions))
            if last_x is not None:
                dx=x-last_x; dy=last_f-f
                denom=dx @ inverse @ dy
                if abs(denom)>1e-20:
                    inverse += np.outer(dx-inverse@dy,dx@inverse)/denom
            last_x=x.copy(); last_f=f.copy()
            step=inverse@f
            scale=min(1.,.15/max(np.linalg.norm(step.reshape(-1,3),axis=1).max(),1e-12))
            x += scale*step
        else:
            raise RuntimeError('Exact V/W joint iteration failed at time %g; no results exported.'%end)
    result.build()
    return result
