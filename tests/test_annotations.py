#===============================================================================
#
#  Reaction and species annotations, and their OMEX Metadata 1.2 export (bg2cellml/annotations.py)
#
#===============================================================================

import asyncio
from pathlib import Path
import sys
import zipfile

import lxml.etree as etree
import pyoxigraph as ox
import pytest

from bg2cellml.annotations import model_annotations, Participant
from bg2cellml.bondgraph.framework import get_framework

#===============================================================================

MODELS = Path(__file__).parent / 'models'
REFERENCES = Path(__file__).parent / 'references'      # the metadata of the test models, which have no annotations

CELLML_URI = 'file:///OUT/model.cellml'
LOCAL = 'file:///OUT/model.ttl#'

BQBIOL = 'http://biomodels.net/biology-qualifiers/'
SEMSIM = 'http://bime.uw.edu/semsim/'
OPB = 'https://identifiers.org/opb:'
XSD_DOUBLE = 'http://www.w3.org/2001/XMLSchema#double'
CMETA_ID = '{http://www.cellml.org/metadata/2.0#}id'

ID = 'https://identifiers.org/'
CYTOSOL = f'{ID}GO:0005829'

#===============================================================================

def model_of(name: str):
    async def load():
        framework = await get_framework()
        source = (MODELS / f'{name}.ttl').resolve()
        model = framework.make_bondgraph_model(source.as_uri(), source.read_text())
        assert not model.has_issues, [issue.reason for issue in model.issues]
        return model
    return asyncio.run(load())

def metadata_of(name: str) -> str:
    cellml_model = model_of(name).make_cellml_model()
    return asyncio.run(cellml_model.metadata(CELLML_URI))

def triples(turtle: str) -> set[tuple[str, str, str]]:
    def term(t):
        return f'"{t.value}"^^{t.datatype.value}' if isinstance(t, ox.Literal) else t.value
    return {(term(t.subject), term(t.predicate), term(t.object)) for t in ox.parse(turtle, format=ox.RdfFormat.TURTLE)}

def objects(graph: set, subject: str, predicate: str) -> set[str]:
    return {o for s, p, o in graph if s == subject and p == predicate}

def fragment(uri: str) -> str:
    return uri.rsplit('#', 1)[-1]

#===============================================================================

def test_reaction_participants_sides_and_coefficients():
    '''Participants are found directly and through 1- and 0-junctions; the bond at the reaction gives the side,
    bondCount the coefficient; a 0-junction doesn't lead on into another reaction's 1-junction.'''
    annotations = model_annotations(model_of('annotated_reactions'))
    reactions = {fragment(r.element): r for r in annotations.reactions}
    assert set(reactions) == {'Re_HK', 'Re_AK', 'Re_XY'}
    participants = lambda name: [p._replace(element=fragment(p.element)) for p in reactions[name].participants]
    assert participants('Re_HK') == [
        Participant('Se_ATP', 'ATP', 'c', 'left', 1),
        Participant('q_Glc', 'Glc', 'c', 'left', 1),
        Participant('Se_ADP', 'ADP', 'c', 'right', 1),
        Participant('Se_H', 'H', 'c', 'right', 1),
        Participant('q_G6P', 'G6P', 'c', 'right', 1),
    ]
    assert participants('Re_AK') == [
        Participant('Se_ADP_2', 'ADP', 'c', 'left', 2),
        Participant('Se_AMP', 'AMP', 'c', 'right', 1),
        Participant('Se_ATP_2', 'ATP', 'c', 'right', 1),
    ]
    assert participants('Re_XY') == [Participant('q_X', 'X', 'c', 'left', 1), Participant('q_Y', 'Y', 'c', 'right', 1)]

def test_reaction_and_species_annotations():
    annotations = model_annotations(model_of('annotated_reactions'))
    hk = next(r for r in annotations.reactions if r.element.endswith('#Re_HK'))
    assert (hk.rhea, hk.rhea_parts, hk.sbo, hk.mediators) == (
        [f'{ID}rhea:17825'], [], [f'{ID}SBO:0000176'], [f'{ID}uniprot:P19367'])
    assert hk.annotated
    assert not any(r.annotated for r in annotations.reactions if not r.element.endswith('#Re_HK'))
    # copies of a chemostat are one entity
    adp = annotations.entities[('ADP', 'c')]
    assert (adp.chebi, adp.go) == ([f'{ID}CHEBI:456216'], [CYTOSOL])
    assert [fragment(e) for e in adp.elements] == ['Se_ADP', 'Se_ADP_2']
    assert annotations.entities[('G6P', 'c')].chebi == [f'{ID}CHEBI:61548']
    assert not annotations.entities[('X', 'c')].annotated

def test_lumped_reaction_parts():
    '''bqbiol:hasPart lists the Rhea reactions a lumped reaction combines.'''
    path = MODELS / 'annotated_reactions.ttl'
    source = path.read_text() + f'<#Re_XY> <{BQBIOL}hasPart> <{ID}rhea:10000>, <{ID}rhea:10004> .\n'
    async def load():
        framework = await get_framework()
        return framework.make_bondgraph_model(path.resolve().as_uri(), source)
    lumped = model_annotations(asyncio.run(load()))
    re_xy = next(r for r in lumped.reactions if r.element.endswith('#Re_XY'))
    assert re_xy.rhea_parts == [f'{ID}rhea:10000', f'{ID}rhea:10004'] and re_xy.rhea == [] and re_xy.annotated

#===============================================================================

def test_the_composite_annotations_exported():
    graph = triples(metadata_of('annotated_reactions'))
    model = 'file:///OUT/model.cellml#'
    process = f'{LOCAL}process_Re_HK'
    # the flux variable is a chemical molar flow rate of the process
    assert objects(graph, f'{model}main-v_Re_HK', f'{BQBIOL}isVersionOf') == {f'{OPB}OPB_00592'}
    assert objects(graph, f'{model}main-v_Re_HK', f'{BQBIOL}isPropertyOf') == {process}
    assert objects(graph, process, f'{BQBIOL}isVersionOf') == {f'{ID}rhea:17825', f'{ID}SBO:0000176'}
    # its sources and sinks, with their multipliers and entities
    def participants(predicate):
        return {(fragment(node), fragment(next(iter(objects(graph, node, f'{SEMSIM}hasPhysicalEntityReference')))),
                 next(iter(objects(graph, node, f'{SEMSIM}hasMultiplier'))))
                for node in objects(graph, process, predicate)}
    one = f'"1"^^{XSD_DOUBLE}'
    assert participants(f'{SEMSIM}hasSourceParticipant') == {
        ('source_Re_HK_1', 'entity_ATP_c', one), ('source_Re_HK_2', 'entity_Glc_c', one)}
    assert participants(f'{SEMSIM}hasSinkParticipant') == {
        ('sink_Re_HK_1', 'entity_ADP_c', one), ('sink_Re_HK_2', 'entity_H_c', one), ('sink_Re_HK_3', 'entity_G6P_c', one)}
    # the enzyme
    mediator = f'{LOCAL}mediator_Re_HK'
    assert objects(graph, process, f'{SEMSIM}hasMediatorParticipant') == {mediator}
    assert objects(graph, mediator, f'{SEMSIM}hasPhysicalEntityReference') == {f'{ID}uniprot:P19367'}
    # entities, and a store's amount as the molar amount of its entity
    assert objects(graph, f'{LOCAL}entity_ATP_c', f'{BQBIOL}is') == {f'{ID}CHEBI:30616'}
    assert objects(graph, f'{LOCAL}entity_ATP_c', f'{BQBIOL}isPartOf') == {CYTOSOL}
    assert objects(graph, f'{model}main-q_Glc_Glc_c', f'{BQBIOL}isVersionOf') == {f'{OPB}OPB_00425'}
    assert objects(graph, f'{model}main-q_Glc_Glc_c', f'{BQBIOL}isPropertyOf') == {f'{LOCAL}entity_Glc_c'}
    # a reaction not itself annotated, of annotated species: bondCount 2 is the multiplier
    source = next(iter(objects(graph, f'{LOCAL}process_Re_AK', f'{SEMSIM}hasSourceParticipant')))
    assert objects(graph, source, f'{SEMSIM}hasMultiplier') == {f'"2"^^{XSD_DOUBLE}'}
    assert objects(graph, f'{LOCAL}process_Re_AK', f'{BQBIOL}isVersionOf') == set()
    # nothing for the reaction and species without annotations
    assert not {s for s, _, _ in graph if s in (f'{LOCAL}process_Re_XY', f'{LOCAL}entity_X_c', f'{LOCAL}entity_Y_c')}
    assert objects(graph, f'{model}main-v_Re_XY', f'{BQBIOL}isPropertyOf') == set()
    assert objects(graph, f'{model}main-q_X_X_c', f'{BQBIOL}isVersionOf') == set()
    # the links to the bond graph elements stay
    assert (f'{model}main-v_Re_HK', 'http://biomodels.net/model-qualifiers/isDerivedFrom',
            (MODELS / 'annotated_reactions.ttl').resolve().as_uri() + '#Re_HK') in graph

def test_annotations_name_the_renamed_variables_by_their_ids():
    '''A naming rule renames the reactions' flows (J_<name>): the annotation refers to them by their cmeta:id,
    which the rename keeps.'''
    cellml_model = model_of('annotated_reactions').make_cellml_model()
    names = {v.get(CMETA_ID): v.get('name')
             for v in etree.fromstring(cellml_model.to_xml().encode()).iter('{http://www.cellml.org/cellml/1.1#}variable')}
    graph = triples(asyncio.run(cellml_model.metadata(CELLML_URI)))
    flows = {fragment(s) for s, p, o in graph if p == f'{BQBIOL}isVersionOf' and o == f'{OPB}OPB_00592'}
    assert {names[id] for id in flows} == {'J_Re_HK', 'J_Re_AK'}

def test_a_model_without_annotations_has_the_same_metadata():
    for name in ('chemostat_copies', 'chemostat_named'):
        metadata = metadata_of(name).replace(MODELS.resolve().as_uri(), 'file:///MODELS')
        reference = (REFERENCES / f'{name}.metadata.ttl').read_text().replace('file:///MODEL.cellml', CELLML_URI)
        assert triples(metadata) == triples(reference)
        assert 'bqbiol' not in metadata and 'local:' not in metadata

#===============================================================================

def test_bgtool_writes_the_annotations(tmp_path, monkeypatch):
    '''--annotate writes them next to the CellML, the OMEX archive holds them with local: resources in the
    annotation file's own namespace.'''
    pytest.importorskip('structlog')
    pytest.importorskip('libopencor')
    monkeypatch.syspath_prepend(str(Path(__file__).parent.parent / 'cli'))
    import bgtool.__main__ as bgtool
    # the CellML's own check runs it with libOpenCOR, which is tested elsewhere
    monkeypatch.setattr(bgtool, 'valid_cellml', lambda cellml: True)
    bg2cellml = bgtool.bg2cellml
    # a name of its own, so the source and the annotation (<stem>.ttl) are different files in the archive
    source = tmp_path / 'annotated_reactions.bgrdf.ttl'
    source.write_text((MODELS / 'annotated_reactions.ttl').read_text())
    source = str(source)
    (tmp_path / 'out').mkdir()
    asyncio.run(bg2cellml(source, str(tmp_path / 'out'), annotate=True, bgrdf=True))
    graph = triples((tmp_path / 'out' / 'annotated_reactions.ttl').read_text())
    local = (tmp_path / 'out' / 'annotated_reactions.ttl').resolve().as_uri() + '#'
    cellml = (tmp_path / 'out' / 'annotated_reactions.cellml').resolve().as_uri() + '#'
    assert objects(graph, f'{cellml}main-v_Re_HK', f'{BQBIOL}isPropertyOf') == {f'{local}process_Re_HK'}

    asyncio.run(bg2cellml(source, str(tmp_path), omex='model.omex', bgrdf=True))
    with zipfile.ZipFile(tmp_path / 'model.omex') as omex:
        archived = omex.read('annotated_reactions.ttl').decode()
    assert '@prefix local: <annotated_reactions.ttl#>' in archived
    graph = triples(f'@base <http://omex-library.org/model.omex/> .\n{archived}')
    assert objects(graph, 'http://omex-library.org/model.omex/annotated_reactions.cellml#main-v_Re_HK',
                   f'{BQBIOL}isPropertyOf') == {'http://omex-library.org/model.omex/annotated_reactions.ttl#process_Re_HK'}

#===============================================================================
#===============================================================================
