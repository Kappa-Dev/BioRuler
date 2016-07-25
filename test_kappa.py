from bioruler.library.kappa_translator import KappaExporter
from regraph.library.data_structures import (TypedDiGraph, TypedHomomorphism)
from regraph.library.utils import plot_graph

def make_action_graph(G, hom):
    res = TypedDiGraph()

    nodes_to_add = [(G.node[n].type_, hom[n], G.node[n].attrs_) for n in G.nodes()]

    edges_to_add = [(G.node[n1].type_, G.node[n2].type_, G.get_edge(n1, n2)) \
                             for n1, n2 in G.edges()]

    for n, t, attrs in nodes_to_add:
        if n in res.nodes():
            res.add_node_attrs(n, attrs)
        else:
            res.add_node(n, t, attrs)
    for n1, n2, attrs in edges_to_add:
        if (n1, n2) in res.edges():
            res.add_edge_attrs(n1, n2, attrs)
        else:
            res.add_edge(n1, n2, attrs)

    return res

graph = TypedDiGraph()
graph.add_nodes_from([
    ('A1', 'A'),
    ('B1', 'B'),

    ('a1_1', 'a_1'),
    ('a1_2', 'a_2'),
    ('b1_1', 'b_1'),
    ('b1_2', 'b_2'),

    ('a1_1_s', 'a_1_s', {'val' : '0'}),

    ('notBND', 'notBND'),
    ('s_BND1', 's_BND1'),
    ('s_BND2', 's_BND2'),

    ('is_BND', 'is_BND'),
    ('s_BND3', 's_BND3'),
    ('s_BND4', 's_BND4'),

    ('MOD1', 'MOD1', {'fun' : {('0', '1')}}),
    ('t_MOD1', 't_MOD1'),
])

graph.add_edges_from([
    ('a1_1', 'A1'), ('a1_2', 'A1'),
    ('a1_1_s', 'a1_1'),
    ('b1_1', 'B1'), ('b1_2', 'B1'),
    ('s_BND1', 'notBND'), ('s_BND2', 'notBND'),
    ('a1_1', 's_BND1'), ('b1_1', 's_BND2'),
    ('s_BND3', 'is_BND'), ('s_BND4', 'is_BND'),
    ('a1_2', 's_BND3'), ('b1_2', 's_BND4'),
    ('t_MOD1', 'MOD1'),
    ('t_MOD1', 'a1_1_s'),
])


action_graph = make_action_graph(
    graph,
    {
        'A1' : 'agent',
        'B1' : 'agent',

        'a1_1' : 'site',
        'a1_2' : 'site',
        'b1_1' : 'site',
        'b1_2' : 'site',

        'a1_1_s' : 'state',

        'notBND' : 'not_BND',
        'is_BND' : 'is_BND',

        's_BND1' : 's_BND',
        's_BND2' : 's_BND',
        's_BND3' : 's_BND',
        's_BND4' : 's_BND',

        'MOD1' : 'MOD',
        't_MOD1' : 't_MOD',
    }
)

action_graph.add_nodes_from([
    ('C', 'agent'),
    ('c_1', 'site'),
    ('BND2', 'BND'),
    ('s_BND21', 's_BND'),
    ('s_BND22', 's_BND'),
    ('BND3', 'BND'),
    ('s_BND31', 's_BND'),
    ('s_BND32', 's_BND'),
])

action_graph.add_edges_from([
    ('c_1', 'C'),
    ('s_BND21', 'BND2'), ('s_BND22', 'BND2'),
    ('s_BND31', 'BND3'), ('s_BND32', 'BND3'),
    ('c_1', 's_BND21'), ('a_1', 's_BND22'),
    ('b_1', 's_BND31'), ('a_1', 's_BND32'),
])

graph.metamodel_ = action_graph
graph.hom = TypedHomomorphism.canonic(graph, action_graph)

print("\n".join(KappaExporter.compile_nugget(graph)))
