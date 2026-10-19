#===============================================================================
#
#  A model's variables named its own way in the CellML (bg2cellml/variable_names.py)
#
#===============================================================================

import asyncio
import json
import math
from pathlib import Path
import re

import pytest

from bg2cellml.bondgraph.framework import get_framework
from bg2cellml.module_export import export_module
from bg2cellml.variable_names import valid_name, variable_table

from .cellml_eval import cellml_from_bgrdf, evaluate

#===============================================================================

MODELS = Path(__file__).parent / 'models'
RT = 8.31446261815324 * 310.0       # the chemical templates' RT

def model_of(name: str):
    async def load():
        framework = await get_framework()
        source = (MODELS / f'{name}.ttl').resolve()
        model = framework.make_bondgraph_model(source.as_uri(), source.read_text())
        assert not model.has_issues, [issue.reason for issue in model.issues]
        return model
    return asyncio.run(load())

def cellml_names(cellml: str) -> set[str]:
    return set(re.findall(r'<variable name="([^"]+)"', cellml))

#===============================================================================

def test_rules_and_explicit_names_rename_the_cellml():
    '''A rule names a variable of a class of element, a rule without a class every element's, an explicit name
    wins over both, and copies of a chemostat named alike get _2.'''
    cellml = cellml_from_bgrdf(MODELS / 'chemostat_named.ttl')
    names = cellml_names(cellml)
    assert {'n_B_c', 'mu_B_c', 'mu_A_c', 'mu_A_c_2', 'J_A_c', 'J_A_c_2', 'v_Re_1', 'v_second',
            'k_Re_1', 'k_Re_2'} <= names
    # without a rule, a variable keeps its bond graph name
    assert 'u_f_Re_1' in names and 'u_r_Re_2' in names
    assert not {'q_B_B_c', 'u_q_B_B_c', 'kappa_Re_1', 'v_Re_2', 'u_Se_A_c', 'v_Se_A_c'} & names
    # a rule for another domain doesn't apply; the rule for the class and domain beats the class's
    assert 'Kd_B' in names and 'dndt_B' in names
    assert not {'never_q_B', 'not_this_B', 'K_q_B_B_c', 'v_q_B_B_c'} & names
    # every use is renamed too: the CellML evaluates as the bond graph's does
    values, rates = evaluate(cellml)
    reference, reference_rates = evaluate(cellml_from_bgrdf(MODELS / 'chemostat_copies.ttl'))
    assert rates['n_B_c'] == pytest.approx(reference_rates['q_B'])
    assert values['v_Re_1'] == pytest.approx(reference['v_Re_1'])
    assert values['v_second'] == pytest.approx(reference['v_Re_2'])
    assert values['mu_A_c'] == values['mu_A_c_2'] == pytest.approx(1000.0)
    # the later copy, by identifier, gets _2
    table = {r['element'].rsplit('#', 1)[-1] + '.' + r['variable']: r['cellml'] for r in variable_table(model_of('chemostat_named'))}
    assert (table['Se_A_1.u'], table['Se_A_2.u']) == ('mu_A_c', 'mu_A_c_2')

def test_a_model_without_naming_keeps_its_names():
    assert cellml_names(cellml_from_bgrdf(MODELS / 'chemostat_copies.ttl')) >= {'q_B', 'u_q_B', 'v_Re_1', 'kappa_Re_2'}

def test_the_table_lists_every_element_variable():
    rows = variable_table(model_of('chemostat_named'))
    by_name = {r['cellml']: r for r in rows}
    assert by_name['n_B_c'] == {'element': by_name['n_B_c']['element'], 'class': 'QuantityStore', 'variable': 'q',
                                'bondgraph': 'q_B_B_c', 'cellml': 'n_B_c',
                                'species': 'B', 'location': 'c', 'name': 'q_B'}
    assert by_name['v_second']['bondgraph'] == 'v_Re_2'
    assert len(rows) == len({r['cellml'] for r in rows})       # all distinct

@pytest.mark.parametrize('name, var, valid', [
    ('mu_H+_c', 'u', 'mu_H_c'),             # a character a name can't hold
    ('3PG_c', 'q', 'q_3PG_c'),              # a digit first
    ('__n__G6P__', 'q', 'n_G6P'),           # repeated and surrounding _
    ('', 'v', 'v'),
])
def test_names_are_made_valid(name, var, valid):
    assert valid_name(name, var) == valid

def test_the_module_export_uses_the_names(tmp_path):
    '''The chemostat's port takes the potential source rule's names; the copy that had the name is renamed.'''
    model = model_of('chemostat_named')
    cellml = model.make_cellml_model().to_xml()
    export_module(model, cellml, tmp_path, 'named', 'bondgraph_v01')
    (entry,) = json.loads((tmp_path / 'named_bondgraph_v01_modules_config.json').read_text())
    assert entry['general_ports'] == [{'port_type': 'chemostat_A_c', 'variables': ['mu_A_c', 'J_A_c']}]
    exported = (tmp_path / 'named_bondgraph_v01_modules.cellml').read_text()
    assert {'mu_A_c', 'J_A_c', 'mu_A_c_copy', 'J_A_c_copy', 'mu_A_c_2', 'J_A_c_2', 'n_B_c', 'n_B_c_init'} \
        <= cellml_names(exported)

#===============================================================================
#===============================================================================
