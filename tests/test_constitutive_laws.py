#===============================================================================
#
#  CellDL and bondgraph tools
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#===============================================================================

"""The CellML bg2cellml makes for elements with the saturating and Bernoulli constitutive relations, a chemical
potential source (a chemostat) and a bond with a bondCount, evaluated at the model's initial values against the
relations worked by hand:

    pytest tests/test_constitutive_laws.py
"""

import math
from pathlib import Path

import pytest

from .cellml_eval import cellml_from_bgrdf, evaluate

MODELS = Path(__file__).parent / 'models'
RT = 2578.73058

#===============================================================================

def test_saturating_reaction_from_a_chemostat_into_a_saturating_store():
    '''2 A -> B: the bond's bondCount 2 doubles A's potential at the reaction; the reaction is mass action over
    1 + (exp(u_f/RT) + exp(u_r/RT))/sigma; B's potential is RT ln(K q k_sat/(k_sat + q)); B gains the reaction's
    flow and the chemostat gives twice it.'''
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'saturating_reaction.ttl'))
    u_a, kappa, sigma, q, K, k_sat = 1000.0, 0.5, 3.0, 0.2, 2.0, 0.5
    u_b = RT*math.log(K*q*k_sat/(k_sat + q))
    pf, pr = math.exp(2*u_a/RT), math.exp(u_b/RT)
    v = kappa*(pf - pr)/(1 + (pf + pr)/sigma)
    assert values['u_f_Re'] == pytest.approx(2*u_a)
    assert values['u_q_B'] == pytest.approx(u_b)
    assert values['v_Re'] == pytest.approx(v)
    assert rates['q_B'] == pytest.approx(v)
    assert abs(values['v_Se_A']) == pytest.approx(2*v)

def test_elements_labelled_alike_keep_their_own_variables():
    '''Two copies of a chemostat with the same symbol, species and location (as a diagram draws one beside each
    reaction) are two sources: each has its own variables, and B gains both reactions' flows.'''
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'chemostat_copies.ttl'))
    u_a, q, K = 1000.0, 0.2, 2.0
    drive = math.exp(u_a/RT) - K*q
    assert values['v_Re_1'] == pytest.approx(0.5*drive)
    assert values['v_Re_2'] == pytest.approx(0.25*drive)
    assert rates['q_B'] == pytest.approx(0.75*drive)
    sources = [name for name in values if name.startswith('v_Se')]
    assert len(sources) == 2

def test_electrogenic_reaction_through_transformers():
    '''A -> B moving one charge into a compartment at phi = Q/C: through two Transformer nodes (between the chemical
    and electrical domains, so not gyrators) the forward potential sees -F phi/2 and the reverse +F phi/2, the reaction
    flows kappa (exp((u_A - F phi/2)/RT) - exp((u_B + F phi/2)/RT)), and the capacitor gains F v.'''
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'electrogenic_reaction.ttl'))
    F = 96485.332
    u_a, kappa, q, K, Q, C = 1000.0, 0.5, 0.2, 2.0, -0.015, 0.1
    phi = Q/C
    v = kappa*(math.exp((u_a - F*phi/2)/RT) - math.exp((RT*math.log(K*q) + F*phi/2)/RT))
    assert values['v_Re'] == pytest.approx(v, rel=1e-6)
    assert rates['q_B'] == pytest.approx(v, rel=1e-6)
    assert rates['q_C_m'] == pytest.approx(F*v, rel=1e-6)

def test_bernoulli_resistor():
    '''A Bernoulli resistor between two hydraulic stores carries the flow whose pressure drop R v + B v|v| is the
    stores' pressure difference.'''
    values, rates = evaluate(cellml_from_bgrdf(MODELS / 'bernoulli_resistor.ttl'))
    R, B = 4.0, 1.5
    drop = values['u_q_1'] - values['u_q_2']
    v = values['v_R']
    assert R*v + B*v*abs(v) == pytest.approx(drop)
    assert rates['q_1'] == pytest.approx(-v) and rates['q_2'] == pytest.approx(v)

#===============================================================================
#===============================================================================
