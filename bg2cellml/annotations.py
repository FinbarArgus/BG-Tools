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

"""
A model's reaction and species annotations, and each reaction's participants.

Reactions (``bgf:Reaction`` elements) are annotated with ``bqbiol:isVersionOf`` a Rhea reaction and an SBO term,
``bqbiol:hasPart`` the Rhea reactions a lumped reaction combines, and ``semsim:hasMediatorParticipant`` a UniProt
protein. Species elements (those with ``bgf:hasSpecies``) are annotated with ``bqbiol:is`` a ChEBI term and
``bqbiol:isPartOf`` a GO cellular component. Terms are identifiers.org IRIs, ``https://identifiers.org/<prefix>:<id>``.

A reaction's participants are the species elements bonded to it, directly or through junctions: through 1-junctions
(the reaction's stoichiometry) and 0-junctions (a species' potential), but not from a 0-junction into a 1-junction,
which belongs to another reaction. The bond at the reaction gives the side: into the reaction is the left
(substrates), out of it the right (products). The coefficient is the product of the bonds' ``bgf:bondCount``
(default 1) along the way.
"""

#===============================================================================

from collections import defaultdict
from dataclasses import dataclass, field
import re
from typing import NamedTuple, TYPE_CHECKING

from .bondgraph.framework_support import ONENODE_JUNCTION, REACTION, ZERONODE_JUNCTION

if TYPE_CHECKING:
    from .bondgraph import BondgraphModel

#===============================================================================

BQBIOL = 'http://biomodels.net/biology-qualifiers/'
SEMSIM = 'http://bime.uw.edu/semsim/'

ANNOTATIONS = f"""
    SELECT DISTINCT ?element ?predicate ?term
    WHERE {{
        <%MODEL%> bgf:hasBondElement ?element .
        ?element ?predicate ?term .
        FILTER (?predicate IN (<{BQBIOL}isVersionOf>, <{BQBIOL}hasPart>, <{BQBIOL}is>, <{BQBIOL}isPartOf>,
                               <{SEMSIM}hasMediatorParticipant>))
    }}"""

SPECIES_ELEMENTS = """
    SELECT DISTINCT ?element ?species ?location
    WHERE {
        <%MODEL%> bgf:hasBondElement ?element .
        ?element bgf:hasSpecies ?species .
        OPTIONAL { ?element bgf:hasLocation ?location }
    }"""

BONDS = """
    SELECT DISTINCT ?bond ?source ?target ?bondCount
    WHERE {
        <%MODEL%> bgf:hasPowerBond ?bond .
        ?bond bgf:hasSource ?source ;
              bgf:hasTarget ?target .
        OPTIONAL { ?bond bgf:bondCount ?bondCount }
    }"""

#===============================================================================

def term_prefix(term: str) -> str|None:
    """The identifiers.org prefix of a term (``rhea``, ``CHEBI``, ``GO``, ``SBO``, ``uniprot``, ...), else None."""
    if (m := re.match(r'https?://identifiers\.org/([^:/]+):', term)) is not None:
        return m.group(1)
    return None

class Participant(NamedTuple):
    element: str            # the species element's URI
    species: str
    location: str
    side: str               # 'left' (a substrate) or 'right' (a product)
    coefficient: int

@dataclass
class Reaction:
    element: str                                            # the reaction element's URI
    rhea: list[str] = field(default_factory=list)           # bqbiol:isVersionOf a Rhea reaction
    rhea_parts: list[str] = field(default_factory=list)     # bqbiol:hasPart, for a lumped reaction
    sbo: list[str] = field(default_factory=list)            # bqbiol:isVersionOf an SBO term
    mediators: list[str] = field(default_factory=list)      # semsim:hasMediatorParticipant
    participants: list[Participant] = field(default_factory=list)

    @property
    def annotated(self) -> bool:
        return bool(self.rhea or self.rhea_parts or self.sbo or self.mediators)

@dataclass
class SpeciesEntity:
    species: str
    location: str
    chebi: list[str] = field(default_factory=list)          # bqbiol:is, from any of its elements
    go: list[str] = field(default_factory=list)             # bqbiol:isPartOf, from any of its elements
    elements: list[str] = field(default_factory=list)       # the species elements with this species and location

    @property
    def annotated(self) -> bool:
        return bool(self.chebi or self.go)

@dataclass
class ModelAnnotations:
    reactions: list[Reaction]                               # every reaction element, in URI order
    entities: dict[tuple[str, str], SpeciesEntity]          # by (species, location)

    def entity(self, participant: Participant) -> SpeciesEntity:
        return self.entities[(participant.species, participant.location)]

#===============================================================================

def model_annotations(model: 'BondgraphModel') -> ModelAnnotations:
    """The model's reactions with their annotations and participants, and its species entities with theirs."""
    uri = model.uri.value
    graph = model.rdf_graph
    terms: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for row in graph.query(ANNOTATIONS.replace('%MODEL%', uri)):
        terms[row['element'].value].append((row['predicate'].value, row['term'].value))   # pyright: ignore[reportOptionalMemberAccess]

    entities: dict[tuple[str, str], SpeciesEntity] = {}
    species_of: dict[str, tuple[str, str]] = {}
    for row in sorted(graph.query(SPECIES_ELEMENTS.replace('%MODEL%', uri)), key=lambda r: r['element'].value):  # pyright: ignore[reportOptionalMemberAccess]
        element = row['element'].value                                                  # pyright: ignore[reportOptionalMemberAccess]
        location = row['location'].value if row.get('location') is not None else ''    # pyright: ignore[reportOptionalMemberAccess]
        key = (row['species'].value, location)                                          # pyright: ignore[reportOptionalMemberAccess]
        species_of[element] = key
        entity = entities.setdefault(key, SpeciesEntity(*key))
        entity.elements.append(element)
        for predicate, term in sorted(terms.get(element, [])):
            prefix = term_prefix(term)
            if predicate == f'{BQBIOL}is' and prefix == 'CHEBI' and term not in entity.chebi:
                entity.chebi.append(term)
            elif predicate == f'{BQBIOL}isPartOf' and prefix == 'GO' and term not in entity.go:
                entity.go.append(term)

    # bonds by node: (bond count, other node, outwards from this node)
    bonds: dict[str, list[tuple[int, str, bool]]] = defaultdict(list)
    for row in graph.query(BONDS.replace('%MODEL%', uri)):
        count = int(row['bondCount'].value) if row.get('bondCount') is not None else 1  # pyright: ignore[reportOptionalMemberAccess]
        source, target = row['source'].value, row['target'].value                       # pyright: ignore[reportOptionalMemberAccess]
        bonds[source].append((count, target, True))
        bonds[target].append((count, source, False))
    junction_types = {junction.uri.value: junction.type for junction in model.junctions}

    def participants(reaction: str) -> list[Participant]:
        found: list[Participant] = []
        for count, node, outwards in bonds.get(reaction, []):
            side = 'right' if outwards else 'left'
            # (node, coefficient so far, the reaction and junctions on the way)
            stack = [(node, count, {reaction})]
            while stack:
                node, coefficient, visited = stack.pop()
                if node in species_of:
                    species, location = species_of[node]
                    found.append(Participant(node, species, location, side, coefficient))
                elif (junction_type := junction_types.get(node)) in (ONENODE_JUNCTION, ZERONODE_JUNCTION):
                    for next_count, next_node, _ in bonds.get(node, []):
                        if next_node in visited:
                            continue
                        if junction_type == ZERONODE_JUNCTION and junction_types.get(next_node) == ONENODE_JUNCTION:
                            continue
                        stack.append((next_node, coefficient*next_count, visited | {node}))
        return sorted(found, key=lambda p: (p.side, p.element))

    reactions: list[Reaction] = []
    for element in sorted(model.elements, key=lambda e: e.uri.value):
        if element.element_class != REACTION:
            continue
        reaction = Reaction(element.uri.value)
        for predicate, term in sorted(terms.get(reaction.element, [])):
            prefix = term_prefix(term)
            if predicate == f'{BQBIOL}isVersionOf' and prefix == 'rhea':
                reaction.rhea.append(term)
            elif predicate == f'{BQBIOL}isVersionOf' and prefix == 'SBO':
                reaction.sbo.append(term)
            elif predicate == f'{BQBIOL}hasPart' and prefix == 'rhea':
                reaction.rhea_parts.append(term)
            elif predicate == f'{SEMSIM}hasMediatorParticipant':
                reaction.mediators.append(term)
        reaction.participants = participants(reaction.element)
        reactions.append(reaction)
    return ModelAnnotations(reactions, entities)

#===============================================================================
#===============================================================================
