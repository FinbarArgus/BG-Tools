#===============================================================================
#
#  The variable mapping of a model's CellML, in libcuflynx's format
#
#===============================================================================

"""
The variable mapping of the CellML that bg2cellml writes, as a CSV in libcuflynx's format
(``libcuflynx.reporting.variable_mapping``, read by ``cuflynx-methods-latex --mapping`` and refreshed by
``cuflynx-variable-mapping``) with one more column:

* ``variable_name``: the canonical ``component/variable`` name (a trailing ``_module`` dropped from the component);
* ``latex``: the symbol in the methods' equations. A variable keeping its bond graph name is written in bond graph
  notation (``q^{\\mathrm{G6P}}_{c}``); one the model names its own way (bgf:NamingRule, bgf:VariableName) by
  libcuflynx's default rule for its name (``c_{\\mathrm{G6P}}``);
* ``kind``: time, state, parameter, constant or algebraic;
* ``component``: the CellML component;
* ``units``: the variable's CellML units;
* ``bondgraph``: the variable's bond graph name (``q_G6P_c``), empty for a variable that isn't an element's.

libcuflynx needs only the first two columns and keeps the others it doesn't know.
"""

import csv
from pathlib import Path
import re
from typing import TYPE_CHECKING

import lxml.etree as etree

from .variable_names import CELLML_NS, MATHML_NS, variable_table

if TYPE_CHECKING:
    from .bondgraph import BondgraphModel

#===============================================================================

COLUMNS = ('variable_name', 'latex', 'kind', 'component', 'units')
BONDGRAPH_COLUMN = 'bondgraph'
MODULE_SUFFIX = '_module'
TIME = 't'

# libcuflynx.reporting.latex_names: the default rule for a name's symbol
GREEK = ('alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu',
         'nu', 'xi', 'pi', 'rho', 'sigma', 'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega')
GREEK_UPPER = ('Gamma', 'Delta', 'Theta', 'Lambda', 'Xi', 'Pi', 'Sigma', 'Upsilon', 'Phi', 'Psi', 'Omega')
_SPECIAL = {'\\': r'\textbackslash{}', '{': r'\{', '}': r'\}', '_': r'\_', '%': r'\%', '$': r'\$', '#': r'\#',
            '&': r'\&', '^': r'\^{}', '~': r'\~{}'}

def _escape(text: str) -> str:
    return ''.join(_SPECIAL.get(ch, ch) for ch in text)

def _greek(token: str) -> str|None:
    if token in GREEK or token in GREEK_UPPER:
        return '\\' + token
    if token.lower() in GREEK and token[0].isupper():
        return '\\' + token.lower()
    return None

def _token(token: str) -> str:
    return _greek(token) or (token if len(token) == 1 else r'\mathrm{' + _escape(token) + '}')

def symbol(name: str) -> str:
    """libcuflynx's default symbol for ``name``: the first part of the name, the rest as its subscripts."""
    tokens = [t for t in name.split('_') if t]
    if not tokens:
        return r'\mathrm{' + _escape(name) + '}'
    head, rest = tokens[0], tokens[1:]
    subs = []
    if (match := re.fullmatch(r'([A-Za-z])(\d+)', head)):
        head_tex = match.group(1)
        subs.append(match.group(2))
    else:
        head_tex = _token(head)
    subs += [_token(t) for t in rest]
    return head_tex + ('_{' + ','.join(subs) + '}' if subs else '')

def bondgraph_symbol(variable: str, species: str, location: str, element: str) -> str:
    """A variable in bond graph notation: a species' as ``q^{species}_{location}``, another element's with the
    element's identifier as a subscript (``v_{\\mathrm{Re\\_HK}}``, ``q_{B}`` for the store q_B)."""
    parts = [p for p in variable.split('_') if p]
    head = _token(parts[0]) if parts else _token(variable)
    subs = [_token(p) for p in parts[1:]]
    if species:
        if location:
            subs.append(_token(location))
        return f'{head}^{{{_token(species)}}}' + ('_{' + ','.join(subs) + '}' if subs else '')
    # the element's identifier, without a leading repeat of the variable (q_B's q is q_{B})
    if parts and element.startswith(parts[0] + '_') and len(element) > len(parts[0]) + 1:
        element = element[len(parts[0]) + 1:]
    subs.append(_token(element))
    return head + '_{' + ','.join(subs) + '}'

#===============================================================================

def _states(component: etree._Element) -> set[str]:
    states = set()
    for apply in component.iter(f'{{{MATHML_NS}}}apply'):
        children = list(apply)
        if children and etree.QName(children[0]).localname == 'diff':
            for child in children[1:]:
                if etree.QName(child).localname == 'ci' and child.text:
                    states.add(child.text.strip())
    return states

def _is_number(text: str|None) -> bool:
    try:
        float(text or '')
        return True
    except ValueError:
        return False

def variable_mapping(model: 'BondgraphModel', cellml: str, parameters: set[str]|None=None) -> list[dict[str, str]]:
    """
    The mapping rows (COLUMNS and ``bondgraph``) of each variable of ``cellml``, the CellML bg2cellml wrote for
    ``model``: time, then each component's states, constants (parameters) and the rest, as libcuflynx orders them.
    ``parameters`` names variables that are a module's parameters (their kind is then parameter).
    """
    parameters = parameters or set()
    table = variable_table(model)
    by_cellml = {row['cellml']: row for row in table}
    root = etree.fromstring(cellml.encode())
    order = {'time': 0, 'state': 1, 'parameter': 2, 'constant': 2, 'algebraic': 3}
    rows = []
    for component in root.iter(f'{{{CELLML_NS}}}component'):
        name = component.get('name', '')
        stripped = name[:-len(MODULE_SUFFIX)] if name.endswith(MODULE_SUFFIX) and len(name) > len(MODULE_SUFFIX) else name
        states = _states(component)
        found = []
        for variable in component.findall(f'{{{CELLML_NS}}}variable'):
            var = variable.get('name', '')
            if var == TIME:
                kind = 'time'
            elif var in states:
                kind = 'state'
            elif var in parameters:
                kind = 'parameter'
            elif _is_number(variable.get('initial_value')):
                kind = 'constant'
            else:
                kind = 'algebraic'
            element = by_cellml.get(var)
            if kind == 'time':
                latex = 't'
            elif element is not None and element['cellml'] == element['bondgraph']:
                latex = bondgraph_symbol(element['variable'], element['species'], element['location'],
                                         element['name'])
            else:
                latex = symbol(var)
            found.append({'variable_name': f'{stripped}/{var}', 'latex': latex, 'kind': kind, 'component': name,
                          'units': variable.get('units', ''),
                          BONDGRAPH_COLUMN: element['bondgraph'] if element is not None else ''})
        rows.extend(sorted(found, key=lambda r: order[r['kind']]))
    return sorted(rows, key=lambda r: r['kind'] != 'time')

def write_variable_mapping(rows: list[dict[str, str]], path: str|Path) -> Path:
    """Write the rows as a CSV of libcuflynx's columns and ``bondgraph``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as fp:
        writer = csv.DictWriter(fp, fieldnames=[*COLUMNS, BONDGRAPH_COLUMN], extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    return path

def mapping_path(cellml_path: str|Path) -> Path:
    """``<dir>/<stem>_variable_mapping.csv`` beside a CellML file, as libcuflynx's default_mapping_path."""
    path = Path(cellml_path)
    stem = path.stem[:-len('_flat')] if path.stem.endswith('_flat') else path.stem
    return path.with_name(f'{stem}_variable_mapping.csv')

#===============================================================================
#===============================================================================
