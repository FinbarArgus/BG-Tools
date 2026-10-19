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

"""Evaluates generated CellML at its initial values without a CellML runtime: every variable with an initial value is
known, and an equation with one unknown left is solved for it, until all are known. Returns the variables' values and
the states' rates."""

import asyncio
from pathlib import Path

import lxml.etree as etree
import sympy

from bg2cellml.bondgraph.framework import get_framework
from bg2cellml.mathml import MATHML_NS, sympy_from_mathml

CELLML_NS = '{http://www.cellml.org/cellml/1.1#}'

#===============================================================================

def cellml_from_bgrdf(path: str|Path) -> str:
#============================================
    async def convert():
        framework = await get_framework()
        assert not framework.has_issues, framework.issues
        source = Path(path).resolve()
        model = framework.make_bondgraph_model(source.as_uri(), source.read_text())
        assert not model.has_issues, [issue.reason for issue in model.issues]
        return model.make_cellml_model().to_xml()
    return asyncio.run(convert())

def evaluate(cellml: str) -> tuple[dict[str, float], dict[str, float]]:
#======================================================================
    root = etree.fromstring(cellml.encode() if isinstance(cellml, str) else cellml)
    known: dict[sympy.Symbol, float] = {}
    for var in root.iter(f'{CELLML_NS}variable'):
        if var.get('initial_value') is not None:
            known[sympy.Symbol(var.get('name'))] = float(var.get('initial_value'))
    rates: dict[str, float] = {}
    pending = []
    for math in root.iter(MATHML_NS.math):
        for eq in math.findall(MATHML_NS.apply):
            _, lhs, rhs = eq.getchildren()
            left, right = sympy_from_mathml(lhs), sympy_from_mathml(rhs)
            if isinstance(left, sympy.Derivative):
                pending.append(('rate', str(left.expr), right))
            else:
                pending.append(('eq', left, right))
    while pending:
        progress = False
        for item in list(pending):
            kind, left, right = item
            if kind == 'rate':
                if right.free_symbols <= set(known):
                    rates[left] = float(right.subs(known))
                    pending.remove(item)
                    progress = True
                continue
            expr = (left - right).subs(known)
            unknown = expr.free_symbols - {sympy.Symbol('t')}
            if not unknown:
                pending.remove(item)
            elif len(unknown) == 1:
                (x,) = unknown
                solution = sympy.solve(expr, x)
                assert len(solution) == 1, (item, solution)
                known[x] = float(solution[0])
                pending.remove(item)
                progress = True
        assert progress, f'cannot solve {pending}'
    return {str(k): v for k, v in known.items()}, rates

#===============================================================================
#===============================================================================
