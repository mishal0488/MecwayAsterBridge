# -*- coding: utf-8 -*-
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Linear releases and exact nonlinear spherical / zero-twist V/W joints."""
import math
from bridge_utils import CONFIG_BRIDGE, get_active_analysis, log


def prepare(md, axis_function):
    if md.get('releases_prepared') or not md.get('beam_releases'):
        return True
    try:
        _prepare(md, axis_function)
        md['releases_prepared'] = True
        return True
    except ValueError as err:
        md['preparation_error'] = 'Beam end releases: ' + str(err)
        log('ERROR: Beam end releases: ' + str(err))
        return False


def _prepare(md, axis_function):
    import bridge_custom
    import bridge_xfem
    analysis = bridge_custom.apply_options(
        get_active_analysis(), md.get('aster_custom', {}).get('options', {}))
    nonlinear = CONFIG_BRIDGE['active_analysis'] == 'nonlin_static_3d'
    if CONFIG_BRIDGE['active_analysis'] not in ('lin_static_3d', 'modal_3d', 'nonlin_static_3d'):
        raise ValueError('supported in linear static, modal and restricted nonlinear static 3D only.')
    if nonlinear and (analysis.get('deformation') != 'GROT_GDEP' or md.get('beam_offsets')):
        raise ValueError('nonlinear releases require GROT_GDEP and no offsets.')
    if CONFIG_BRIDGE.get('result_route', 'liml') != 'liml':
        raise ValueError('select the LIML results route.')
    if (md.get('shell_upgrades') or
            CONFIG_BRIDGE.get('elem_hex8') == 'LINE_QUAD' or
            md.get('contacts') or md.get('bonded') or md.get('couplings') or
            md.get('cyclic') or analysis.get('adapt_refine') or
            analysis.get('modal_complex') or bridge_xfem.is_requested(md)):
        raise ValueError('mesh upgrades, contact/couplings, cyclic, complex modes, adaptive refinement and XFEM are not supported together with releases.')
    elements = dict((e[0], e) for e in md['mesh_elements'])
    selected = {}
    allowed = ('POU_D_T_GD',) if nonlinear else ('POU_D_E', 'POU_D_T')
    for item in md['beam_releases']:
        faces = md.get('face_selection_faces', {}).get(item['selection'])
        if not faces:
            raise ValueError('release selection %r must contain beam end faces.' % item['selection'])
        for eid, face in faces:
            e = elements.get(eid)
            if (not e or e[1] != 'line2' or face not in (2, 3) or
                    md['elset_modelisations'].get(e[3]) not in allowed):
                raise ValueError('release on element %s must be an end face of a line2 %s beam.' % (eid, '/'.join(allowed)))
            selected.setdefault((eid, face - 2), set()).update(item['axes'])
    selected = dict((k, v) for k, v in selected.items() if v)
    if nonlinear and any(axes not in (set(('u','v','w')), set(('v','w'))) for axes in selected.values()):
        raise ValueError('geometric nonlinear analysis supports full U/V/W releases and combined V/W releases only; other selective releases need different finite-rotation constraints.')
    custom_vw = nonlinear and any(axes == set(('v','w')) for axes in selected.values())
    if custom_vw:
        if any(md['elset_modelisations'].get(e[3]) != 'POU_D_T_GD'
               for e in md['mesh_elements']):
            raise ValueError('custom nonlinear V/W joints currently require a beam-only POU_D_T_GD model.')
        custom = md.get('aster_custom', {})
        if (analysis.get('crit_stab') or analysis.get('imperfection_amplitude') or
                md.get('thermal_stress') is not None or
                any(custom.get(k) for k in ('code','excit','post','failure','fatigue','welding','phases','metallurgy'))):
            raise ValueError('custom nonlinear V/W joints support elastic static mechanical loading only, without stability extraction, automatic imperfections, thermal histories or custom analysis code.')
    # Retained reference-node rotations must belong to a real, unreleased
    # beam end. Otherwise unused rotational DOFs need separate elimination.
    unreleased = set()
    for eid, shape, nodes, group in md['mesh_elements']:
        if md['elset_modelisations'].get(group, '').startswith('POU_'):
            for end, nid in enumerate(nodes):
                if (eid, end) not in selected:
                    unreleased.add(nid)
    for eid, end in selected:
        if elements[eid][2][end] not in unreleased:
            raise ValueError('released node %s must also connect to an unreleased beam end; isolated/all-released junctions are not supported yet.' % elements[eid][2][end])
    for bc in md['bcs']:
        faces = md.get('face_selection_faces', {}).get(bc.get('selection'), [])
        if any((eid, face-2) in selected for eid, face in faces):
            raise ValueError('use node selections for loads/supports on released beam ends.')
    coords = dict(md['mesh_nodes'])
    new = max(coords) + 1
    aliases, equations, changed, vw_joints = {}, [], [], []
    for eid, shape, nodes, group in md['mesh_elements']:
        ends = list(nodes)
        if any((eid, end) in selected for end in range(len(nodes))):
            p, q = [coords[n] for n in nodes]
            length = math.sqrt(sum((q[i]-p[i])**2 for i in range(3)))
            if not length:
                raise ValueError('zero-length beam %s.' % eid)
            u = tuple((q[i]-p[i])/length for i in range(3))
            v = axis_function(p, q, md['beam_orientations'].get(eid))
            # VECT_Y is projected by Code_Aster too (notably for Mecway's
            # nearly-vertical default). Use the same orthonormal frame.
            dot = sum(u[i]*v[i] for i in range(3))
            v = tuple(v[i]-dot*u[i] for i in range(3))
            size = math.sqrt(sum(c*c for c in v))
            if size < 1.e-12:
                raise ValueError('undefined local axes for beam %s.' % eid)
            v = tuple(c/size for c in v)
            w = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
            for end, old in enumerate(nodes):
                released = selected.get((eid, end))
                if not released:
                    continue
                coords[new] = coords[old]
                aliases[new] = old
                ends[end] = new
                for dof in ('DX', 'DY', 'DZ'):
                    equations.append({'terms': [(new, dof, 1.), (old, dof, -1.)], 'constant': 0.})
                is_vw = nonlinear and released == set(('v','w'))
                if is_vw:
                    vw_joints.append({'reference': old, 'end': new, 'axis': u})
                for name, axis in zip(('u', 'v', 'w'), (u, v, w)):
                    if is_vw:
                        continue  # Nonlinear quaternion constraint, never a fixed-axis tie.
                    if name not in released:
                        terms = []
                        for dof, c in zip(('DRX', 'DRY', 'DRZ'), axis):
                            if abs(c) > 1.e-14:
                                terms.extend([(new, dof, c), (old, dof, -c)])
                        equations.append({'terms': terms, 'constant': 0.})
                new += 1
        changed.append((eid, shape, ends, group))
    md['mesh_nodes'] = coords
    md['mesh_elements'] = changed
    # The existing internal-node result mapping also handles coincident
    # released ends: element identity is retained and reactions are summed.
    md['offset_node_aliases'] = aliases
    md['release_equations'] = equations
    if vw_joints:
        md['nonlinear_vw_joints'] = vw_joints
        log('INFO: Custom exact nonlinear V/W joints: translations coincide, two swing freedoms, zero relative swing-twist torsion. Only converged joint states are exported; check load-step refinement.')
    elif nonlinear:
        log('INFO: Exact nonlinear spherical joints: coincident translations tied; all relative rotations free.')
    log('INFO: Beam end releases: %s independent beam ends; translations and unreleased local rotations tied exactly.' % len(aliases))
