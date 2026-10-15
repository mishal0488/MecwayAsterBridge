# -*- coding: utf-8 -*-
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Small-rotation beam offsets, compatible with IronPython 2.7.

Physical centrelines carry stiffness and mass. Reference nodes carry nodal
loads/supports; exact linear rigid-arm equations connect the two.
"""
import math
import re
from bridge_utils import CONFIG_BRIDGE, get_active_analysis, log, parse_and_convert


def read_vector(item):
    """Reject expressions/unknown units instead of silently treating them as zero."""
    values = []
    for name in ('u', 'v', 'w'):
        children = item.findall(name)
        text = item.findtext(name, '0 m').strip()
        match = re.match(r'^[-+]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][-+]?[0-9]+)?\s*(m|mm|cm|km|in|inch|ft|um)?$', text)
        if not match or len(children) > 1 or (children and children[0].attrib):
            values.append(float('nan'))
        else:
            values.append(parse_and_convert(text, 'length'))
    return tuple(values)


def prepare(md, axis_function):
    if md.get('offsets_prepared') or not md.get('beam_offsets'):
        return True
    try:
        _prepare(md, axis_function)
        md['offsets_prepared'] = True
        return True
    except ValueError as err:
        log('ERROR: Beam offsets: ' + str(err))
        return False


def _prepare(md, axis_function):
    import bridge_custom
    import bridge_xfem
    analysis = bridge_custom.apply_options(
        get_active_analysis(), md.get('aster_custom', {}).get('options', {}))
    name = CONFIG_BRIDGE['active_analysis']
    if name not in ('lin_static_3d', 'modal_3d'):
        raise ValueError('supported only in linear static and modal 3D; large-rotation offsets are not implemented.')
    if CONFIG_BRIDGE.get('result_route', 'liml') != 'liml':
        raise ValueError('select the LIML results route.')
    if CONFIG_BRIDGE.get('elem_hex8') == 'LINE_QUAD' or md.get('shell_upgrades'):
        raise ValueError('mesh order upgrades are not supported with offsets.')
    if md.get('contacts') or md.get('bonded') or md.get('couplings') or md.get('cyclic'):
        raise ValueError('contact, bonded/coupled or cyclic models are not yet supported.')
    if analysis.get('adapt_refine') or bridge_xfem.is_requested(md):
        raise ValueError('adaptive refinement and XFEM are not supported with offsets.')
    if name == 'modal_3d' and (analysis.get('modal_prestress') or analysis.get('modal_complex')):
        raise ValueError('disable modal prestress and complex modes for offset beams.')
    elements = dict((e[0], e) for e in md['mesh_elements'])
    selected = {}
    for item in md['beam_offsets']:
        ids = md['element_selections'].get(item['selection'])
        if ids is None:
            ids = [e[0] for e in elements.values() if e[3] == item['selection']]
        if not ids:
            raise ValueError('offset selection %r is empty or is not an element selection.' % item['selection'])
        for eid in ids:
            if eid in selected:
                raise ValueError('overlapping offset definitions at element %s.' % eid)
            e = elements.get(eid)
            if not e or e[1] != 'line2' or md['elset_modelisations'].get(e[3]) not in ('POU_D_E', 'POU_D_T'):
                raise ValueError('only two-node POU_D_E/POU_D_T beams may have offsets.')
            if not all(isinstance(x, (int, float)) and not math.isnan(x) and not math.isinf(x) for x in item['vector']):
                raise ValueError('offsets must be finite constant lengths.')
            selected[eid] = item['vector']
    for bc in md['bcs']:
        selection = bc.get('selection')
        affected = set(md.get('element_selections', {}).get(selection, []))
        affected.update(eid for eid, fid in md.get('face_selection_faces', {}).get(selection, []))
        affected.update(e[0] for e in elements.values() if e[3] == selection)
        if affected.intersection(selected) and bc.get('type') != 'gravity':
            raise ValueError('use node selections for loads/supports in offset models (selection %r).' % selection)
    coords = dict(md['mesh_nodes'])
    next_node = max(coords) + 1
    aliases, refs, equations, changed = {}, set(), [], []
    for eid, shape, nodes, group in md['mesh_elements']:
        if eid not in selected or not any(selected[eid]):
            changed.append((eid, shape, nodes, group))
            continue
        p, q = [coords[n] for n in nodes]
        length = math.sqrt(sum((q[i]-p[i])**2 for i in range(3)))
        if not length:
            raise ValueError('zero-length beam %s.' % eid)
        u = tuple((q[i]-p[i])/length for i in range(3))
        v = axis_function(p, q, md['beam_orientations'].get(eid))
        w = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
        r = tuple(sum(selected[eid][j]*a[i] for j, a in enumerate((u, v, w))) for i in range(3))
        new_nodes = []
        for old in nodes:
            new = next_node
            next_node += 1
            coords[new] = tuple(coords[old][i]+r[i] for i in range(3))
            new_nodes.append(new)
            aliases[new] = old
            refs.add(old)
            # u_new = u_reference + theta_reference cross r.
            cross = ((0, r[2], -r[1]), (-r[2], 0, r[0]), (r[1], -r[0], 0))
            for i, dof in enumerate(('DX', 'DY', 'DZ')):
                terms = [(new, dof, 1.0), (old, dof, -1.0)]
                terms += [(old, d, -cross[i][j]) for j, d in enumerate(('DRX', 'DRY', 'DRZ')) if cross[i][j]]
                equations.append({'terms': terms, 'constant': 0.0})
            for dof in ('DRX', 'DRY', 'DRZ'):
                equations.append({'terms': [(new, dof, 1.0), (old, dof, -1.0)], 'constant': 0.0})
        changed.append((eid, shape, new_nodes, group))
    md['mesh_nodes'] = coords
    md['mesh_elements'] = changed
    md['offset_node_aliases'] = aliases
    md['offset_reference_nodes'] = sorted(refs)
    md['offset_equations'] = equations
    log('INFO: Beam offsets: %s reference nodes, %s physical nodes; linear rigid arms with no artificial stiffness or mass.' % (len(refs), len(aliases)))
