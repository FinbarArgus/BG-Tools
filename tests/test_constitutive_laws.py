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
