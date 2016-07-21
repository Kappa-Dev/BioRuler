"""."""
from bioruler.library.metamodels import metamodel_kappa

from regraph.library.data_structures import (TypedDiGraph,
                                             Homomorphism,
                                             TypedHomomorphism)
from regraph.library.utils import (keys_by_value,
                                   plot_graph)

import networkx as nx
import subprocess
import os
import json

from copy import deepcopy


class KappaExporter(object):

    meta = metamodel_kappa

    def __init__(self):
        pass

    @staticmethod
    def check_nugget(G, hom=None):

        # check syntax
        if G.metamodel_ is None:
            raise ValueError(
                "You didn't provide the action graph"
            )
        else:
            if hom is None:
                hom = TypedHomomorphism.canonic(G.metamodel_, KappaExporter.meta)
            else:
                if hom.source_ != G.metamodel_:
                    raise ValueError(
                        "Invalid homomorphism, source is not the action graph "+\
                        "of G"
                    )
                elif hom.target_ != KappaExporter.meta:
                    raise ValueError(
                        "Invalid homomorphism, target is not our metamodel"
                    )
                elif set(hom.mapping_.values()) > set(KappaExporter.meta.nodes()):
                    raise ValueError(
                        "Invalid homomorphism, mapping is not consistent"
                    )

        # check semantic

        # BRK in action graph
        BRK_a_g = keys_by_value(hom.mapping_, 'BRK')
        # BRK in nugget
        BRK_g = set()
        for brk in BRK_a_g:
            BRK_g.update(keys_by_value(G.hom.mapping_, brk))

        # BRK should have exactly 2 entering edges to a t_BRK node

        for brk in BRK_g:
            in_edges = G.in_edges(brk)
            if len(in_edges) != 2:
                raise ValueError(
                    "Node %s of type BRK have %s entering edges, expected 2"%
                    (brk, len(in_edges))
                )
            t1 = in_edges[0][0]
            t2 = in_edges[1][0]

            if hom.mapping_[G.hom.mapping_[t1]] != 't_BRK':
                raise ValueError(
                    "Node %s have type %s, expected type t_BRK" %
                    (t1, hom.mapping_[G.hom.mapping_[t1]])
                )
            if hom.mapping_[G.hom.mapping_[t2]] != 't_BRK':
                raise ValueError(
                    "Node %s have type %s, expected type t_BRK" %
                    (t2, hom.mapping_[G.hom.mapping_[t2]])
                )

        # BND in action graph
        BND_a_g = keys_by_value(hom.mapping_, 'BND')
        # BND in nugget
        BND_g = set()
        for bnd in BND_a_g:
            BND_g.update(keys_by_value(G.hom.mapping_, bnd))

        # BND should have exactly 2 entering edges to a s node

        for bnd in BND_g:
            in_edges = G.in_edges(bnd)
            if len(in_edges) != 2:
                raise ValueError(
                    "Node %s of type BND have %s entering edges, expected 2"%
                    (bnd, len(in_edges))
                )
            s1 = in_edges[0][0]
            s2 = in_edges[1][0]

            if hom.mapping_[G.hom.mapping_[s1]] != 's_BND':
                raise ValueError(
                    "Node %s have type %s, expected type s" %
                    (s1, hom.mapping_[G.hom.mapping_[s1]])
                )
            if hom.mapping_[G.hom.mapping_[s2]] != 's_BND':
                raise ValueError(
                    "Node %s have type %s, expected type s" %
                    (s2, hom.mapping_[G.hom.mapping_[s2]])
                )

        # is_BND should have exactly 2 entering edges to a s node

        # isBND in action graph
        isBND_a_g = keys_by_value(hom.mapping_, 'is_BND')
        # isBND in nugget
        isBND_g = set()
        for isbnd in isBND_a_g:
            isBND_g.update(keys_by_value(G.hom.mapping_, isbnd))

        for isbnd in isBND_g:
            in_edges = G.in_edges(isbnd)
            if len(in_edges) > 2:
                raise ValueError(
                    "Node %s of type is_BND have %s entering edges, expected 2 at most"%
                    (isbnd, len(in_edges))
                )
            for s in in_edges:
                s = s[0]

                if hom.mapping_[G.hom.mapping_[s]] != 's_BND':
                    raise ValueError(
                        "Node %s have type %s, expected type s" %
                        (s, hom.mapping_[G.hom.mapping_[s]])
                    )

        # MOD should have exactly 1 entering edge to a t_MOD node

        # MOD in action graph
        MOD_a_g = keys_by_value(hom.mapping_, 'MOD')
        # MOD in nugget
        MOD_g = set()
        for mod in MOD_a_g:
            MOD_g.update(keys_by_value(G.hom.mapping_, mod))

        for mod in MOD_g:
            in_edges = G.in_edges(mod)
            if len(in_edges) != 1:
                raise ValueError(
                    "Node %s of type MOD have %s entering edges, expected 1"%
                    (brk, len(in_edges))
                )
            t = in_edges[0][0]

            if hom.mapping_[G.hom.mapping_[t]] != 't_MOD':
                raise ValueError(
                    "Node %s have type %s, expected type t_MOD" %
                    (t, hom.mapping_[G.hom.mapping_[t]])
                )

        # state should have exactly 1 leaving edges to a site

        # state in action graph
        state_a_g = keys_by_value(hom.mapping_, 'state')
        # state in nugget
        state_g = set()
        for state in state_a_g:
            state_g.update(keys_by_value(G.hom.mapping_, state))

        for state in state_g:
            out_edges = G.out_edges(state)
            t_mod_nodes = []

            if len(out_edges) != 1:
                raise ValueError(
                    "Node %s of type state have %s leaving edges, expected 1" %
                    (state, len(out_edges))
                )

            if hom.mapping_[G.hom.mapping_[out_edges[0][1]]] != 'site':
                raise ValueError(
                    "Node %s have type %s, expected it to be site because of edge "%
                    (out_edges[0][1], hom.mapping_[G.hom.mapping_[out_edges[0][1]]])+\
                    "%s-%s" % out_edges[0]
                )

            for node, state in G.in_edges(state):
                if hom.mapping_[G.hom.mapping_[node]] == 't_MOD':
                    t_mod_nodes.append(node)

            if len(t_mod_nodes) > 1 :
                raise ValueError(
                    "Inconsistent nugget, state %s of node %s can't be modified " %
                    (state, agent) +\
                    " twice"
                )

        # site should have exactly 1 leaving edge to an agent

        # site in action graph
        site_a_g = keys_by_value(hom.mapping_, 'site')
        # site in nugget
        site_g = set()
        for site in site_a_g:
            site_g.update(keys_by_value(G.hom.mapping_, site))

        for site in site_g:
            out_edges = G.out_edges(site)
            have_agent = False
            s_nodes = []
            t_brk_nodes = []
            for n1, n2 in out_edges:
                if hom.mapping_[G.hom.mapping_[n2]] == 'agent':
                    if have_agent == True:
                        have_agent = None
                    else:
                        have_agent = True
                        agent = n2
                elif hom.mapping_[G.hom.mapping_[n2]] == 's_BND':
                    s_nodes.append(n2)

            if have_agent is None:
                raise ValueError(
                    "Node %s of type site is linked to more than one agent" % site
                )
            if not have_agent:
                raise ValueError(
                    "Node %s of type site isn't linked with any agent" % site
                )
            if len(s_nodes) > 1 :
                raise ValueError(
                    "Inconsistent nugget, site %s of node %s can't bind twice, " %
                    (site, agent) +\
                    "or maybe you have a redundant is_BND node."
                )

            for node, site in G.in_edges(site):
                if hom.mapping_[G.hom.mapping_[node]] == 't_BRK':
                    t_brk_nodes.append(node)

            if len(t_brk_nodes) > 1 :
                raise ValueError(
                    "Inconsistent nugget, site %s of node %s can't unbind twice, " %
                    (site, agent)
                )

        # syn/deg should have at least 1 entering edge from a t_sd or a s_sd

        # syn/deg in action graph
        sd_a_g = keys_by_value(hom.mapping_, 'SYN/DEG')
        # syn/deg in nugget
        sd_g = set()
        for sd in sd_a_g:
            sd_g.update(keys_by_value(G.hom.mapping_, sd))

        for sd in sd_g:
            if len(G.in_edges(sd)) < 1:
                raise ValueError(
                    "Node %s of type SYN/DEG doesn't have any source nor target"%
                    sd
                )

    @staticmethod
    def normalize_nugget(graph, hom=None):

        def merge_states(G, states):
            vals = []
            for s in states:
                if vals == []:
                    vals = [str(v) for v in G.node[s].attrs_['val']]
                else:
                    new_val = []
                    for old_val in vals:
                        for val in G.node[s].attrs_['val']:
                            new_val.append(old_val+str(val))
                    if new_val != []:
                        vals = new_val
            name = G.merge_nodes(states)
            G.node[name].attrs_['val'] = set(vals)

            return name

        # check syntax
        if graph.metamodel_ is None:
            raise ValueError(
                "You didn't provide the action graph"
            )
        else:
            if hom is None:
                hom = TypedHomomorphism.canonic(graph.metamodel_, KappaExporter.meta)
            else:
                if hom.source_ != graph.metamodel_:
                    raise ValueError(
                        "Invalid homomorphism, source is not the action graph "+\
                        "of G"
                    )
                elif hom.target_ != KappaExporter.meta:
                    raise ValueError(
                        "Invalid homomorphism, target is not our metamodel"
                    )
                elif set(hom.mapping_.values()) > set(KappaExporter.meta.nodes()):
                    raise ValueError(
                        "Invalid homomorphism, mapping is not consistent"
                    )


        res = deepcopy(graph)

        # Merge states of action graph

        # site in action graph
        site_g = keys_by_value(hom.mapping_, 'site')

        for site in site_g:
            if site in graph.hom.mapping_.values():
                states = []
                for state, site in res.metamodel_.in_edges(site):
                    if hom.mapping_[state] == 'state':
                        states.append(state)
                if len(states) > 1:
                    s_pred = []
                    for s in states :
                        s_pred+= keys_by_value(res.hom.mapping_, s)

                    name = merge_states(res.metamodel_, states)

                    for el in s_pred:
                        res.hom.mapping_[el] = name
                        res.node[el].type_ = name

                    for s in states :
                        del hom.mapping_[s]
                    hom.mapping_[name] = 'state'

        # Merge nodes of nugget

        # site in action graph
        site_a_g = keys_by_value(hom.mapping_, 'site')
        # site in nugget
        site_g = []
        for site in site_a_g:
            site_ = keys_by_value(res.hom.mapping_, site)
            site_g += site_

        for site in site_g:
            states = []
            for state, site in res.in_edges(site):
                if hom.mapping_[res.hom.mapping_[state]] == 'state':
                    states.append(state)

            if len(states) > 1:
                name = merge_states(res, states)

                site_a_g = res.hom.mapping_[site]
                state_a_g = res.metamodel_.in_edges(site_a_g)[0][0]
                for s in states :
                    del res.hom.mapping_[s]
                res.hom.mapping_[name] = state_a_g

        # BND in action_graph
        bnd_a_g = keys_by_value(hom.mapping_, 'BND')
        # BND in nugget
        bnd_g = []
        for bnd in bnd_a_g:
            bnd_ = keys_by_value(res.hom.mapping_, bnd)
            bnd_g += bnd_

        site_bnd = []
        for i in range(len(bnd_g)):
            bnd = bnd_g[i]

            s1 = res.in_edges(bnd)[0][0]
            s2 = res.in_edges(bnd)[1][0]

            aux = []
            for site1, s1 in res.in_edges(s1):
                for site2, s2 in res.in_edges(s2):
                    aux.append((site1, site2))

            site_bnd.append(aux)


        # BRK in action_graph
        brk_a_g = keys_by_value(hom.mapping_, 'BRK')
        # BRK in nugget
        brk_g = []
        for brk in brk_a_g:
            brk_ = keys_by_value(res.hom.mapping_, brk)
            brk_g += brk_

        site_brk = []
        for i in range(len(brk_g)):
            brk = brk_g[i]

            t1 = res.in_edges(brk)[0][0]
            t2 = res.in_edges(brk)[1][0]

            aux = []
            for site1, t1 in res.in_edges(t1):
                for site2, t2 in res.in_edges(t2):
                    aux.append((site1, site2))

            site_brk.append(aux)

        # isBND in action graph
        isbnd_a_g = keys_by_value(hom.mapping_, 'isBND')
        # isBND in nugget
        isbnd_g = []
        for isbnd in isbnd_a_g:
            isbnd_ = keys_by_value(res.hom.mapping_, isbnd)
            isbnd_g += isbnd_

        site_isbnd = []
        for i in range(len(isbnd_g)):
            isbnd = isbnd_g[i]

            s1 = res.in_edges(isbnd)[0][0]
            if len(res.in_edges(isbnd)) > 1:
                s2 = res.in_edges(isbnd)[1][0]

            aux = []
            if len(res.in_edges(isbnd)) > 1:
                for site1, s1 in res.in_edges(s1):
                    for site2, s2 in res.in_edges(s2):
                        aux.append((site1, site2))
            else:
                for site1, s1 in res.in_edges(s1):
                    aux.append(site1)

            sites_isbnd.append(aux)

        elements_to_remove = []

        for el in site_brk:
            for s1, s2 in el:
                for i in range(len(site_isbnd)):
                    el = site_isbnd[i]
                    if type(el) == tuple:
                        if s1 == el[0] or s2 == el[0] or\
                           s1 == el[1] or s2 == el[1]:
                           if i not in elements_to_remove:
                               elements_to_remove.append(i)
                    else:
                        if s1 == el or s2 == el:
                            if i not in elements_to_remove:
                                elements_to_remove.append(i)

        nodes_to_remove = []
        for i in elements_to_remove:
            isbnd = isbnd_g[i]
            nodes_to_remove.append(isbnd)
            for s, isbnd in res.in_edges(isbnd):
                nodes_to_remove.append(s)

        for n in nodes_to_remove:
            res.remove_node(n)

        return res

    @staticmethod
    def compile(G, hom=None, agent_sites={}):
        """ G : nugget
            hom : action_graph -> meta_model
            The nugget G have to follow some conditions :
                - It's syntaxically correct : It projects to our meta_model
                - It's semantically correct :
                    - We have two sources for each BND, is_BND, two targets for each
                      BRK and one target for each MOD, each state is linked to one
                      site exactly, each site is linked to one agent exactly
                    - We don't have redundancies (e.g is_BND and BRK linked to the
                      same sites)
                    - We don't have inconsistancies (e.g BND A.s->B.s and A.s->C.s)
                - It's fully expanded : we don't have more than one edge for each
                  source or target node
        """


        # check syntax
        if G.metamodel_ is None:
            raise ValueError(
                "You didn't provide the action graph"
            )
        else:
            if hom is None:
                hom = TypedHomomorphism.canonic(G.metamodel_, KappaExporter.meta)
            else:
                if hom.source_ != G.metamodel_:
                    raise ValueError(
                        "Invalid homomorphism, source is not the action graph "+\
                        "of G"
                    )
                elif hom.target_ != KappaExporter.meta:
                    raise ValueError(
                        "Invalid homomorphism, target is not our metamodel"
                    )
                elif set(hom.mapping_.values()) > set(KappaExporter.meta.nodes()):
                    raise ValueError(
                        "Invalid homomorphism, mapping is not consistent"
                    )


        KappaExporter.check_nugget(G, hom)
        G = KappaExporter.normalize_nugget(G, hom)

        rules=[({}, {})]

        def do_SD(a_sources, a_targets, LHS, RHS):

            for a in a_sources:

                if a not in LHS.keys():
                    LHS[a] = {}

                for site, a in G.in_edges(a):
                    if hom.mapping_[G.hom.mapping_[site]] == 'site':
                        state = None
                        for st, site in G.in_edges(site):
                            if hom.mapping_[G.hom.mapping_[st]] == 'state':
                                if G.node[st].attrs_['val'] != set():
                                    state = deepcopy(G.node[st].attrs_['val']).pop()

                        # updating LHS

                        if site not in LHS[a].keys():
                            LHS[a][site] = [state, None]

            for a in a_targets:

                if a not in RHS.keys():
                    RHS[a] = {}

                for site, a in G.in_edges(a):
                    if hom.mapping_[G.hom.mapping_[site]] == 'site':
                        state = None
                        for st, site in G.in_edges(site):
                            if hom.mapping_[G.hom.mapping_[st]] == 'state':
                                if G.node[st].attrs_['val'] != set():
                                    state = deepcopy(G.node[st].attrs_['val']).pop()

                        # updating RHS

                        if site not in RHS[a].keys():
                            RHS[a][site] = [state, 'isFree']

        def do_BRK(site1, site2, LHS, RHS):
                # target agents
                for e in G.out_edges(site1):
                    if hom.mapping_[G.hom.mapping_[e[1]]] == 'agent':
                        a1 = e[1]
                for e in G.out_edges(site2):
                    if hom.mapping_[G.hom.mapping_[e[1]]] == 'agent':
                        a2 = e[1]

                if a1 not in LHS.keys():
                    LHS[a1] =  {}
                if a2 not in LHS.keys():
                    LHS[a2] = {}

                if a1 not in RHS.keys():
                    RHS[a1] = {}
                if a2 not in RHS.keys():
                    RHS[a2] = {}

                # checking for specified states

                state1 = None
                state2 = None
                for node, site1 in G.in_edges(site1):
                    if hom.mapping_[G.hom.mapping_[node]] == 'state':
                        if G.node[node].attrs_['val'] != set():
                            state1 = deepcopy(G.node[node].attrs_['val'])
                            state1 = state1.pop()
                for node, site2 in G.in_edges(site2):
                    if hom.mapping_[G.hom.mapping_[node]] == 'state':
                        if G.node[node].attrs_['val'] != set():
                            state2 = deepcopy(G.node[node].attrs_['val'])
                            state2 = state2.pop()

                # updating LHS of rule

                LHS[a1][site1] = [state1,
                                           "%s.%s_%s.%s" %
                                           (str(a1), str(site1),
                                            str(a2), str(site2))]

                LHS[a2][site2] = [state2,
                                           "%s.%s_%s.%s" %
                                           (str(a1), str(site1),
                                            str(a2), str(site2))]

                # updating RHS of rule

                RHS[a1][site1] = [state1, None]

                RHS[a2][site2] = [state2, None]

        def do_BND(site1, site2, LHS, RHS):

            # source agents
            for e in G.out_edges(site1):
                if hom.mapping_[G.hom.mapping_[e[1]]] == 'agent':
                    a1 = e[1]
            for e in G.out_edges(site2):
                if hom.mapping_[G.hom.mapping_[e[1]]] == 'agent':
                    a2 = e[1]

            if a1 not in LHS.keys() and a1 not in RHS.keys():
                LHS[a1] = {}
            if a1 not in RHS.keys():
                RHS[a1] = {}
            if a2 not in LHS.keys() and a2 not in RHS.keys():
                LHS[a2] = {}
            if a2 not in RHS.keys():
                RHS[a2] = {}

            # checking for specified states

            state1 = None
            state2 = None
            for node, site1 in G.in_edges(site1):
                if hom.mapping_[G.hom.mapping_[node]] == 'state':
                    if G.node[node].attrs_['val'] != set():
                        state1 = deepcopy(G.node[node].attrs_['val'])
                        state1 = state1.pop()
            for node, site2 in G.in_edges(site2):
                if hom.mapping_[G.hom.mapping_[node]] == 'state':
                    if G.node[node].attrs_['val'] != set():
                        state2 = deepcopy(G.node[node].attrs_['val'])
                        state2 = state2.pop()

            # updating LHS of rule

            if a1 in LHS.keys():
                if site1 not in LHS[a1].keys():
                    LHS[a1][site1] = [state1, None]
            if a2 in LHS.keys():
                if site2 not in LHS[a2].keys():
                    LHS[a2][site2] = [state2, None]

            # updating RHS of rule

            RHS[a1][site1] = [state1,
                                       "%s.%s_%s.%s" %
                                       (str(a1), str(site1),
                                        str(a2), str(site2))]

            RHS[a2][site2] = [state2,
                                       "%s.%s_%s.%s" %
                                       (str(a1), str(site1),
                                        str(a2), str(site2))]

        def do_isBND(site1, site2, LHS, RHS):
            # source agents
            for e in G.out_edges(site1):
                if hom.mapping_[G.hom.mapping_[e[1]]] == 'agent':
                    a1 = e[1]
            for e in G.out_edges(site2):
                if hom.mapping_[G.hom.mapping_[e[1]]] == 'agent':
                    a2 = e[1]

            if a1 not in LHS.keys():
                LHS[a1] = {}
            if a2 not in LHS.keys():
                LHS[a2] = {}

            if a1 not in RHS.keys():
                RHS[a1] = {}
            if a2 not in RHS.keys():
                RHS[a2] = {}

            # checking for specified states

            state1 = None
            state2 = None
            for node, site1 in G.in_edges(site1):
                if hom.mapping_[G.hom.mapping_[node]] == 'state':
                    if G.node[node].attrs_['val'] != set():
                        state1 = deepcopy(G.node[node].attrs_['val'])
                        state1 = state1.pop()
            for node, site2 in G.in_edges(site2):
                if hom.mapping_[G.hom.mapping_[node]] == 'state':
                    if G.node[node].attrs_['val'] != set():
                        state2 = deepcopy(G.node[node].attrs_['val'])
                        state2 = state2.pop()

            # updating LHS of rule

            LHS[a1][site1] = [state1,
                                       "%s.%s_%s.%s" %
                                       (str(a1), str(site1),
                                        str(a2), str(site2))]

            LHS[a2][site2] = [state2,
                                       "%s.%s_%s.%s" %
                                       (str(a1), str(site1),
                                        str(a2), str(site2))]

            # updating RHS of rule

            if site1 not in RHS[a1].keys():
                RHS[a1][site1] = [state1,
                                  "%s.%s_%s.%s" %
                                  (str(a1), str(site1),
                                   str(a2), str(site2))]
            if site2 not in RHS[a2].keys():
                RHS[a2][site2] = [state2,
                                  "%s.%s_%s.%s" %
                                  (str(a1), str(site1),
                                   str(a2), str(site2))]

        def do_MOD(state, LHS, RHS):
            # site targeted
            site = G.out_edges(state)[0][1]

            # agent targeted
            for e in G.out_edges(site):
                if hom.mapping_[G.hom.mapping_[e[1]]] == 'agent':
                    a = e[1]

            if a not in LHS.keys():
                LHS[a] = {}

            if a not in RHS.keys():
                RHS[a] = {}

            if G.node[state].attrs_['val'] == set():
                old_state = None
            else:
                old_state = deepcopy(G.node[state].attrs_['val'])
                old_state = old_state.pop()
            new_state = G.node[mod].attrs_['fun']
            if type(new_state) != set:
                new_state = set([new_state])
            for el in new_state:
                if type(el) == tuple:
                    if el[0] == old_state:
                        res = el[1]
                else:
                    res = el
            new_state = res

            # updating LHS of rule

            if site not in LHS[a].keys():
                LHS[a][site] = [old_state, None]

            # updating RHS of rule

            if site in RHS[a].keys():
                RHS[a][site][0] = new_state
            else:
                RHS[a][site] = [new_state, None]

        def do_isFREE(site, LHS, RHS):
            # agent targeted
            for e in G.out_edges(site):
                if hom.mapping_[G.hom.mapping_[e[1]]] == 'agent':
                    a = e[1]

            if a not in LHS.keys():
                LHS[a] = {}

            if a not in RHS.keys():
                RHS[a] = {}

            state = None
            for node, site in G.in_edges(site):
                if hom.mapping_[G.hom.mapping_[node]] == 'state':
                    if G.node[node].attrs_['val'] != set():
                        state = deepcopy(G.node[node].attrs_['val'])
                        state = state.pop()

            # updating LHS of rule

            LHS[a][site] = [state, "isFree"]

            # updating RHS of rule

            if site not in RHS[a].keys():
                RHS[a][site] = [state, "isFree"]

        A_G = G.metamodel_

        # finding agents

        # agent_site : { agent : { site1 : [states],
        #                          site2 : [states],
        #                          ...}
        #                ...}
        #               default = {}

        agents_a_g = keys_by_value(hom.mapping_, 'agent')
        agents = []
        for a in agents_a_g:
            agents += keys_by_value(G.hom.mapping_, a)

        for a in agents:
            if a in G.hom.mapping_.keys():
                a = G.hom.mapping_[a]
                if not a in agent_sites.keys():
                    agent_sites[a] = {}
                for e in A_G.in_edges(a):
                    if hom.mapping_[e[0]] == 'site':
                        s = e[0]
                        agent_sites[a][s] = []

                        for e_site in A_G.in_edges(s):
                            if hom.mapping_[e_site[0]] == 'state':
                                if A_G.node[e_site[0]].attrs_ is not None and\
                                   'val' in A_G.node[e_site[0]].attrs_.keys():
                                    if len(agent_sites[a][s]) == 0 :
                                        agent_sites[a][s] = list(A_G.node[e_site[0]].attrs_['val'])
                                    else:
                                        new_val = []
                                        for old_val in agent_sites[a][s]:
                                            for v in A_G.node[e_site[0]].attrs_['val']:
                                                new_val.append(str(old_val)+str(v))
                                        agent_sites[a][s] = new_val
                                else:
                                    raise ValueError(
                                        ("State %s of site %s of agent %s doesn't have" %
                                        (e_site[0], s, a))+" 'val' attribute"
                                    )


        # making rules

        # LHS : { agent : { site1 : [state, binding],
        #                   site2 : [state, binding],
        #                   ...}
        #         ...}

        # RHS : { agent : { site1 : [state, binding],
        #                   site2 : [state, binding],
        #                   ...}
        #         ...}

        # rules : [(LHS, RHS), ...]
        #         default = [({}, {})]


        # BRK in action graph
        BRK_a_g = keys_by_value(hom.mapping_, 'BRK')
        # BRK in nugget
        BRK_g = set()
        for brk in BRK_a_g:
            BRK_g.update(keys_by_value(G.hom.mapping_, brk))

        for brk in BRK_g:
            # targets of BRK
            t1 = G.in_edges(brk)[0][0]
            t2 = G.in_edges(brk)[1][0]

            aux_rules = []

            for t1, site1 in G.out_edges(t1):
                for t2, site2 in G.out_edges(t2):
                    if hom.mapping_[G.hom.mapping_[site1]] == 'site' and\
                       hom.mapping_[G.hom.mapping_[site2]] == 'site':
                        aux_aux_rules = deepcopy(rules)
                        for LHS, RHS in aux_aux_rules:
                            do_BRK(site1, site2, LHS, RHS)
                        aux_rules += aux_aux_rules

            rules = aux_rules

        # syn/deg in action graph
        sd_a_g = keys_by_value(hom.mapping_, 'SYN/DEG')
        # syn/deg in nugget
        sd_g = set()
        for sd in sd_a_g:
            sd_g.update(keys_by_value(G.hom.mapping_, sd))

        for sd in sd_g:
            sources = []
            targets = []
            for s, sd in G.in_edges(sd):
                if hom.mapping_[G.hom.mapping_[s]] == 's_SD':
                    sources.append(s)
                else:
                    targets.append(s)
            agent_sources = []
            for s in sources:
                aux = []
                if agent_sources == []:
                    for a, s in G.in_edges(s):
                        aux.append([a])
                else:
                    for a, s in G.in_edges(s):
                        aux_aux = agent_sources
                        for i in range(len(aux_aux)):
                            aux_aux[i].append(a)
                    aux += aux_aux
                agent_sources = aux

            agent_targets = []
            for s in targets:
                aux = []
                if agent_targets == []:
                    for s, a in G.out_edges(s):
                        if hom.mapping_[G.hom.mapping_[a]] == 'agent':
                            aux.append([a])
                else:
                    for s, a in G.out_edges(s):
                        if hom.mapping_[G.hom.mapping_[a]] == 'agent':
                            aux_aux = agent_targets
                            for i in range(len(aux_aux)):
                                aux_aux[i].append(a)
                    aux += aux_aux
                agent_targets = aux

            aux_rules = []

            if agent_targets == []:
                for source in agent_sources:
                    aux_aux_rules = deepcopy(rules)
                    for LHS, RHS in aux_aux_rules:
                        do_SD(source, [], LHS, RHS)
                    aux_rules += aux_aux_rules
            elif agent_sources == []:
                for target in agent_targets:
                    aux_aux_rules = deepcopy(rules)
                    for LHS, RHS in aux_aux_rules:
                        do_SD([], target, LHS, RHS)
                    aux_rules += aux_aux_rules
            else:
                for source in agent_sources:
                    for target in agent_targets:
                        aux_aux_rules = deepcopy(rules)
                        for LHS, RHS in aux_aux_rules:
                            do_SD(source, target, LHS, RHS)
                        aux_rules += aux_aux_rules


            rules = aux_rules



        # BND in action graph
        BND_a_g = keys_by_value(hom.mapping_, 'BND')
        # BND in nugget
        BND_g = set()
        for bnd in BND_a_g:
            BND_g.update(keys_by_value(G.hom.mapping_, bnd))
        for bnd in BND_g:
            # sources of BND
            s1 = G.in_edges(bnd)[0][0]
            s2 = G.in_edges(bnd)[1][0]

            aux_rules = []

            for site1, t1 in G.in_edges(s1):
                for site2, t2 in G.in_edges(s2):
                    aux_aux_rules = deepcopy(rules)
                    for LHS, RHS in aux_aux_rules:
                        do_BND(site1, site2, LHS, RHS)
                    aux_rules += aux_aux_rules

            rules = aux_rules

        # is_BND in action graph
        isBND_a_g = keys_by_value(hom.mapping_, 'is_BND')
        # is_BND in nugget
        isBND_g = set()
        for isbnd in isBND_a_g:
            isBND_g.update(keys_by_value(G.hom.mapping_, isbnd))

        for isbnd in isBND_g:
            # sources of is_BND
            s1 = G.in_edges(isbnd)[0][0]
            if len(G.in_edges(isbnd)) > 1:
                s2 = G.in_edges(isbnd)[1][0]

            aux_rules = []

            if len(G.in_edges(isbnd)) > 1:
                for site1, t1 in G.in_edges(s1):
                    for site2, t2 in G.in_edges(s2):
                        aux_aux_rules = deepcopy(rules)
                        for LHS, RHS in aux_aux_rules:
                            do_isBND(site1, site2, LHS, RHS)
                        aux_rules += aux_aux_rules
            else:
                for site1, t1 in G.in_edges(s1):
                    aux_aux_rules = deepcopy(rules)
                    for LHS, RHS in aux_aux_rules:
                        # source agents
                        for e in G.out_edges(site1):
                            if hom.mapping_[G.hom.mapping_[e[1]]] == 'agent':
                                a1 = e[1]

                        if a1 not in LHS.keys():
                            LHS[a1] = {}

                        if a1 not in RHS.keys():
                            RHS[a1] = {}

                        # checking for specified states

                        state1 = None
                        for node, site1 in G.in_edges(site1):
                            if hom.mapping_[G.hom.mapping_[node]] == 'state':
                                if G.node[node].attrs_['val'] != set():
                                    state1 = deepcopy(G.node[node].attrs_['val'])
                                    state1 = state1.pop()

                        # updating LHS of rule

                        LHS[a1][site1] = [state1,
                                          "_"]

                        # updating RHS of rule

                        if site1 not in RHS[a1].keys():
                            RHS[a1][site1] = [state1, "_"]
                    aux_rules += aux_aux_rules




            rules = aux_rules

        # MOD in action graph
        MOD_a_g = keys_by_value(hom.mapping_, 'MOD')
        # MOD in nugget
        MOD_g = set()
        for mod in MOD_a_g:
            MOD_g.update(keys_by_value(G.hom.mapping_, mod))

        for mod in MOD_g:
            # target of MOD
            t = G.in_edges(mod)[0][0]

            aux_rules = []

            for t, state in G.out_edges(t):
                if hom.mapping_[G.hom.mapping_[state]] == 'state':
                    aux_aux_rules = deepcopy(rules)
                    for LHS, RHS in aux_aux_rules:
                        do_MOD(state, LHS, RHS)
                    aux_rules += aux_aux_rules

            rules = aux_rules

        # is_free in action graph
        isfree_a_g = keys_by_value(hom.mapping_, 'is_FREE')
        # is_free in nugget
        isfree_g = set()
        for isfree in isfree_a_g:
            isfree_g.update(keys_by_value(G.hom.mapping_, isfree))

        for isfree in isfree_g:
            # target of is_FREE
            t = G.in_edges(isfree)[0][0]

            aux_rules = []

            for t, site in G.out_edges(t):
                if hom.mapping_[G.hom.mapping_[site]] == 'site':
                    aux_aux_rules = deepcopy(rules)
                    for LHS, RHS in aux_aux_rules:
                        do_isFREE(site, LHS, RHS)
                    aux_rules += aux_aux_rules

            rules = aux_rules


        context = {}
        for a in G.nodes():
            for LHS, RHS in rules:
                if a not in LHS.keys() and\
                   a not in RHS.keys() and\
                   hom.mapping_[G.hom.mapping_[a]] == "agent":

                    context[a] = {}

                    for site, a in G.in_edges(a):
                        if hom.mapping_[G.hom.mapping_[site]] == 'site':
                            state = None
                            for st, site in G.in_edges(site):
                                if hom.mapping_[G.hom.mapping_[st]] == 'state':
                                    if G.node[st].attrs_['val'] != set():
                                        state = deepcopy(G.node[st].attrs_['val']).pop()

                            # updating context

                            context[a][site] = [state, None]



        return agent_sites, rules, context

    @staticmethod
    def generate_agent_decl(agent_sites):
        # generating agent declarations

        agent_decl = ''

        for a in agent_sites.keys():
            agent_decl += '%agent: '+str(a)+'('
            for s in agent_sites[a].keys():
                agent_decl += str(s)
                for state in agent_sites[a][s]:
                    agent_decl += '~'+str(state)
                agent_decl += ','
            if agent_decl[-1] == ',' : agent_decl = agent_decl[:-1]
            agent_decl += ')\n'

        return agent_decl

    @staticmethod
    def generate_rules_decl(G, rules, context={}, count=0):
        # generating rules

        rule_no = 1
        rules_decl = ''

        for LHS, RHS in rules:
            if (LHS != {} or RHS != {}) and\
                LHS != RHS:
                i = 0
                bindings = {'_':'_'}
                rule = ''
                for a in context.keys():
                    rule += str(G.hom.mapping_[a])+'('
                    for site in context[a].keys():
                        rule += str(G.hom.mapping_[site])
                        state, binding = context[a][site]
                        if state is not None:
                            rule += '~'+str(state)
                        if binding is not None:
                            if binding != 'isFree':
                                if binding in bindings.keys():
                                    rule += '!'+bindings[binding]
                                else:
                                    rule += '!'+str(i)
                                    bindings[binding] = str(i)
                                    i += 1
                        else:
                            rule += '?'
                        rule += ', '
                    if rule[-2:] == ', ' : rule = rule[:-2]
                    rule += '), '
                for a in LHS.keys() :
                    rule += str(G.hom.mapping_[a])+'('
                    for site in LHS[a].keys():
                        rule += str(G.hom.mapping_[site])
                        state, binding = LHS[a][site]
                        if state is not None:
                            rule += '~'+str(state)
                        if binding is not None:
                            if binding != 'isFree':
                                if binding in bindings.keys():
                                    rule += '!'+bindings[binding]
                                else:
                                    rule += '!'+str(i)
                                    bindings[binding] = str(i)
                                    i += 1
                        else:
                            rule += '?'
                        rule += ', '
                    if rule[-2:] == ', ' : rule = rule[:-2]
                    rule += '), '

                if rule[-2:] == ', ' : rule = rule[:-2]

                rule+= ' -> '

                for a in context.keys():
                    rule += str(G.hom.mapping_[a])+'('
                    for site in context[a].keys():
                        rule += str(G.hom.mapping_[site])
                        state, binding = context[a][site]
                        if state is not None:
                            rule += '~'+str(state)
                        if binding is not None:
                            if binding != 'isFree':
                                if binding in bindings.keys():
                                    rule += '!'+bindings[binding]
                                else:
                                    rule += '!'+str(i)
                                    bindings[binding] = str(i)
                                    i += 1
                        else:
                            rule += '?'
                        rule += ', '
                    if rule[-2:] == ', ' : rule = rule[:-2]
                    rule += '), '
                for a in RHS.keys() :
                    rule += str(G.hom.mapping_[a])+'('
                    for site in RHS[a].keys():
                        rule += str(G.hom.mapping_[site])
                        state, binding = RHS[a][site]
                        if state is not None:
                            rule += '~'+str(state)
                        if binding is not None:
                            if binding != 'isFree':
                                if binding in bindings.keys():
                                    rule += '!'+bindings[binding]
                                else:
                                    rule += '!'+str(i)
                                    bindings[binding] = str(i)
                                    i += 1
                        else:
                            rule += '?'
                        rule += ', '
                    if rule[-2:] == ', ' : rule = rule[:-2]
                    rule += '), '
                if rule[-2:] == ', ' : rule = rule[:-2]

                rule += " @ 'NUG_%s_RULE_%s'" % (count, rule_no)
                rule_no += 1
                rules_decl += rule+"\n"

        return rules_decl

    @staticmethod
    def compile_graph(G, hom=None):


        # check syntax
        if G.metamodel_ is None:
            raise ValueError(
                "You didn't provide the action graph"
            )
        else:
            if hom is None:
                hom = TypedHomomorphism.canonic(G.metamodel_, KappaExporter.meta)
            else:
                if hom.source_ != G.metamodel_:
                    raise ValueError(
                        "Invalid homomorphism, source is not the action graph "+\
                        "of G"
                    )
                elif hom.target_ != KappaExporter.meta:
                    raise ValueError(
                        "Invalid homomorphism, target is not our metamodel"
                    )
                elif set(hom.mapping_.values()) > set(KappaExporter.meta.nodes()):
                    raise ValueError(
                        "Invalid homomorphism, mapping is not consistent"
                    )


        undir_G = nx.Graph()
        undir_G.add_nodes_from(G.nodes())
        undir_G.add_edges_from(G.edges())
        con_comp = nx.connected_components(undir_G)
        agent_sites={}
        rules=[]
        for n_nugget in con_comp:
            nugget = G.subgraph(n_nugget)
            agent_sites, rule, context = KappaExporter.compile(nugget, hom, agent_sites)
            rules.append((rule, context))
            agent_decl = KappaExporter.generate_agent_decl(agent_sites)
        rules_decl = ''
        count = 1
        for rule, context in rules:
            rules_decl += KappaExporter.generate_rules_decl(G, rule, context, count)
            count += 1

        return agent_decl, rules_decl

    @staticmethod
    def compile_nugget(nug, hom=None):
        agent_sites, rule, context = KappaExporter.compile(nug, hom)
        agent_decl = KappaExporter.generate_agent_decl(agent_sites)
        rules_decl = KappaExporter.generate_rules_decl(nug, rule, context, 1)

        return agent_decl, rules_decl

    @staticmethod
    def compile_nugget_list(nug_list):
        agent_sites = {}
        rules = ''
        count=1

        for nug in nug_list:
            agent_sites, rule, context = KappaExporter.compile(nug, None, agent_sites)
            rules += KappaExporter.generate_rules_decl(nug, rule, context, count)
            count+=1


        agent_decl = KappaExporter.generate_agent_decl(agent_sites)

        return agent_decl, rules

        return agent_decl, rules_decl


class KappaImporter(object):

    def __init__(self):
        pass

    @staticmethod
    def uncompile_agents(agents):
        action_graph = TypedDiGraph()
        for agent, sites in agents.items():
            action_graph.add_node(agent, "agent")
            for site, states in sites.items():

                action_graph.add_node(agent+"-"+site, "site")
                action_graph.add_edge(agent+"-"+site, agent)

                action_graph.add_node(agent+"-"+site+"-s",
                                      "state",
                                      attrs={'val':set(states)})
                action_graph.add_edge(agent+"-"+site+"-s", agent+"-"+site)

        return action_graph

    @staticmethod
    def uncompile_rule(rule, graph, action_graph, count='0'):
        def gen_name(name, graph=graph):
            i = 0
            new_name = name+"-"+count
            while new_name in graph.nodes():
                i+=1
                new_name = name+"-"+count+"-"+str(i)
            return new_name

        def gen_node(name, name_a_g, link, graph=graph):
            target_name = gen_name(name, graph)
            graph.add_node(target_name, name_a_g)
            graph.add_edge(target_name, link)

            return target_name
        LHS, RHS = rule
        aux = {}
        for k,v in LHS:
            new_name = k
            while new_name in aux.keys():
                new_name += "."+k
            aux[new_name] = {'name':k, 'sites':v}
        LHS = aux

        aux = {}
        for k,v in RHS:
            new_name = k
            while new_name in aux.keys():
                new_name += "."+k
            aux[new_name] = {'name':k, 'sites':v}
        RHS = aux

        agents = []
        a_sites = []
        agents_names = []
        sites_names = []
        el_count = {}

        for a, s in LHS.items():
            agents.append(a)
            a_sites.append(list(s['sites'].items()))

        for a, s in RHS.items():
            if a not in agents:
                agents.append(a)
                a_sites.append(list(s['sites'].items()))
            else:
                y = agents.index(a)
                for site, v in a_sites[y]:
                    if site not in s['sites'].keys():
                        agents.append(a)
                        a_sites.append(list(s['sites'].items()))

        for i in range(len(agents)):
            agent = agents[i]
            sites = a_sites[i]

            if agent in LHS.keys():
                agent_type = LHS[agent]['name']
            else:
                agent_type = RHS[agent]['name']

            agent_name = gen_name(agent_type)
            agents_names.append(agent_name)
            graph.add_node(agent_name, agent_type)

            names = []
            for site, value in sites:
                site_name = gen_name(agent_type+"-"+site)
                graph.add_node(site_name, agent_type+"-"+site)
                names.append(site_name)
                graph.add_edge(site_name, agent_name)

                state = value["state"]
                if state != []:
                    state_name = site_name+"-s"
                    graph.add_node(state_name,
                                   agent_type+"-"+site+"-s",
                                   attrs={'val':state[0]})
                    graph.add_edge(state_name, site_name)
            sites_names.append(names)

        treated_sd = [False for i in agents]

        # SYN/DEG

        creations = []
        deletions = []
        for i in range(len(agents)):
            a1 = agents[i]

            if a1 not in LHS.keys():
                creations.append((a1,agents_names[i]))
                treated_sd[i] = True
            elif a1 not in RHS.keys():
                deletions.append((a1,agents_names[i]))
                treated_sd[i] = True
            else:
                a1_indexes = [i for i, x in enumerate(agents) if x == a1 if not treated_sd[i]]
                if len(a1_indexes) > 1:
                    a1_ind1 = a1_indexes[0]
                    a1_ind2 = a1_indexes[1]
                    if a_sites[a1_ind1] != a_sites[a1_ind2]:
                        if LHS[a1]['sites'] == a_sites[a1_ind1]:
                            creations.append((a1, agents_names[a1_ind1]))
                            deletions.append((a1, agents_names[a1_ind2]))
                        else:
                            creations.append((a1, agents_names[a1_ind2]))
                            deletions.append((a1, agents_names[a1_ind1]))
                        treated_sd[a1_ind1] = True
                        treated_sd[a1_ind2] = True

        if creations+deletions != []:
            syn_name_a_g = gen_name("SYN/DEG", action_graph)
            action_graph.add_node(syn_name_a_g, "SYN/DEG")

            source_name_a_g = gen_name("s_SD", action_graph)
            action_graph.add_node(source_name_a_g, "s_SD")
            action_graph.add_edge(source_name_a_g, syn_name_a_g)

            target_name_a_g = gen_name("t_SD", action_graph)
            action_graph.add_node(target_name_a_g, "t_SD")
            action_graph.add_edge(target_name_a_g, syn_name_a_g)

            syn_name = gen_name("SYN/DEG")
            graph.add_node(syn_name, syn_name_a_g)

        for node in deletions:

            node_type = LHS[node[0]]['name']

            source_name = gen_name("s_SD")
            graph.add_node(source_name, source_name_a_g)
            graph.add_edge(source_name, syn_name)

            action_graph.add_edge(node_type, source_name_a_g)
            graph.add_edge(node[1], source_name)

        for node in creations:

            node_type = RHS[node[0]]['name']

            target_name = gen_name("t_SD")
            graph.add_node(target_name, target_name_a_g)
            graph.add_edge(target_name, syn_name)

            new_node = node[1]

            action_graph.add_edge(target_name_a_g, node_type)
            graph.add_edge(target_name, new_node)

        bounds = []

        for i in range(len(agents)):
            a1 = agents[i]

            if a1 in LHS.keys():
                a1_type = LHS[a1]['name']
            else:
                a1_type = RHS[a1]['name']

            if a1 in LHS.keys() and a1 in RHS.keys():

                for k in range(len(a_sites[i])):
                    # MOD
                    s1 = a_sites[i][k][0]
                    if s1 in RHS[a1]['sites'].keys() and\
                       RHS[a1]['sites'][s1]['state'] != []:
                       if s1 in LHS[a1]['sites'].keys() and\
                          LHS[a1]['sites'][s1]["state"] !=\
                          RHS[a1]['sites'][s1]["state"]:

                              new_val = RHS[a1]['sites'][s1]["state"][0]

                              mod_name_a_g = gen_name("MOD", action_graph)

                              if LHS[a1]['sites'][s1]["state"] != []:
                                  old_val = LHS[a1]['sites'][s1]["state"][0]

                                  action_graph.add_node(mod_name_a_g,
                                                        "MOD",
                                                        attrs={
                                                         'fun': (old_val, new_val)
                                                        })
                              else:
                                  action_graph.add_node(mod_name_a_g,
                                                        "MOD",
                                                        attrs={
                                                         'fun': new_val
                                                        })

                              target_name_a_g = gen_node("t_MOD", "t_MOD", mod_name_a_g, action_graph)

                              mod_name = gen_name("MOD")
                              if LHS[a1]['sites'][s1]["state"] != []:
                                  old_val = LHS[a1]['sites'][s1]["state"][0]

                                  graph.add_node(mod_name,
                                                 mod_name_a_g,
                                                 attrs={
                                                  'fun': (old_val, new_val)
                                                 })
                              else:
                                  graph.add_node(mod_name,
                                                 mod_name_a_g,
                                                 attrs={
                                                  'fun': new_val
                                                 })

                              target_name = gen_node("t_MOD", target_name_a_g, mod_name)

                              action_graph.add_edge(target_name_a_g, a1_type+"-"+s1+"-s")

                              if sites_names[i][k]+"-s" not in graph.nodes():
                                  graph.add_node(sites_names[i][k]+"-s", a1_type+"-"+s1+"-s", {'val':set()})
                                  graph.add_edge(sites_names[i][k]+"-s", sites_names[i][k])

                              graph.add_edge(target_name, sites_names[i][k]+"-s")

                    # is_FREE
                    if s1 in LHS[a1]['sites'].keys() and\
                       LHS[a1]['sites'][s1]['binding'] is None:

                        isfree_name_a_g = gen_name("is_FREE", action_graph)
                        action_graph.add_node(isfree_name_a_g, "is_FREE")

                        target_name_a_g = gen_node("t_FREE", "t_FREE", isfree_name_a_g, action_graph)

                        isfree_name = gen_name("is_FREE")
                        graph.add_node(isfree_name, isfree_name_a_g)

                        target_name = gen_node("t_FREE", target_name_a_g, isfree_name)

                        action_graph.add_edge(target_name_a_g, a1_type+"-"+s1)
                        graph.add_edge(target_name, sites_names[i][k])

                    # is_BND 1 elt
                    if s1 in LHS[a1]['sites'].keys() and\
                       LHS[a1]['sites'][s1]['binding'] == '_':

                        isbnd_name_a_g = gen_name("is_BND", action_graph)
                        action_graph.add_node(isbnd_name_a_g, "is_BND")

                        sbnd1_name_a_g = gen_node("s_BND", "s_BND", isbnd_name_a_g, action_graph)

                        isbnd_name = gen_name("is_BND")
                        graph.add_node(isbnd_name, isbnd_name_a_g)

                        sbnd1_name = gen_node("s_BND", sbnd1_name_a_g, isbnd_name)

                        action_graph.add_edge(a1_type+"-"+s1, sbnd1_name_a_g)
                        graph.add_edge(sites_names[i][k], sbnd1_name)

                for j in range(i+1, len(agents)):
                    a2 = agents[j]
                    if a2 in LHS.keys() and a2 in RHS.keys():

                        if a2 in LHS.keys():
                            a2_type = LHS[a2]['name']
                        else:
                            a2_type = RHS[a2]['name']

                        for k in range(len(a_sites[i])):
                            for l in range(k, len(a_sites[j])):
                                s1 = a_sites[i][k][0]
                                s2 = a_sites[j][l][0]
                                if s1 in LHS[a1]['sites'].keys() and\
                                   s1 in RHS[a1]['sites'].keys() and\
                                   s2 in LHS[a2]['sites'].keys() and\
                                   s2 in RHS[a2]['sites'].keys():

                                    old_bind1 = LHS[a1]['sites'][s1]["binding"]
                                    old_bind2 = LHS[a2]['sites'][s2]["binding"]

                                    new_bind1 = RHS[a1]['sites'][s1]["binding"]
                                    new_bind2 = RHS[a2]['sites'][s2]["binding"]

                                    if old_bind1 == old_bind2 and\
                                       old_bind1 is not None:
                                      if new_bind1 == old_bind1 and\
                                          new_bind2 == old_bind2:

                                          # isBND 2 elts

                                          isbnd_name_a_g = gen_name("is_BND", action_graph)
                                          action_graph.add_node(isbnd_name_a_g, "is_BND")

                                          sbnd1_name_a_g = gen_node("s_BND", "s_BND", isbnd_name_a_g, action_graph)
                                          sbnd2_name_a_g = gen_node("s_BND", "s_BND", isbnd_name_a_g, action_graph)

                                          isbnd_name = gen_name("is_BND")
                                          graph.add_node(isbnd_name, isbnd_name_a_g)

                                          sbnd1_name = gen_node("s_BND", sbnd1_name_a_g, isbnd_name)
                                          sbnd2_name = gen_node("s_BND", sbnd2_name_a_g, isbnd_name)

                                          action_graph.add_edge(a1_type+"-"+s1, sbnd1_name_a_g)
                                          action_graph.add_edge(a2_type+"-"+s2, sbnd2_name_a_g)
                                          graph.add_edge(sites_names[i][k], sbnd1_name)
                                          graph.add_edge(sites_names[j][l], sbnd2_name)

                                      else:

                                          if old_bind1 != "_":

                                              # BRK

                                              brk_name_a_g = gen_name("BRK", action_graph)
                                              action_graph.add_node(brk_name_a_g, "BRK")

                                              tbrk1_name_a_g = gen_node("t_BRK", "t_BRK", brk_name_a_g, action_graph)
                                              tbrk2_name_a_g = gen_node("t_BRK", "t_BRK", brk_name_a_g, action_graph)

                                              brk_name = gen_name("BRK")
                                              graph.add_node(brk_name, brk_name_a_g)

                                              tbrk1_name = gen_node("t_BRK", tbrk1_name_a_g, brk_name)
                                              tbrk2_name = gen_node("t_BRK", tbrk2_name_a_g, brk_name)

                                              action_graph.add_edge(tbrk1_name_a_g, a1_type+"-"+s1)
                                              action_graph.add_edge(tbrk2_name_a_g, a2_type+"-"+s2)
                                              graph.add_edge(tbrk1_name, sites_names[i][k])
                                              graph.add_edge(tbrk2_name, sites_names[j][l])
                                    else:
                                        if new_bind1 == new_bind2:

                                            # BND

                                            bnd_name_a_g = gen_name("BND", action_graph)
                                            action_graph.add_node(bnd_name_a_g, "BND")

                                            sbnd1_name_a_g = gen_node("s_BND", "s_BND", bnd_name_a_g, action_graph)
                                            sbnd2_name_a_g = gen_node("s_BND", "s_BND", bnd_name_a_g, action_graph)

                                            bnd_name = gen_name("BND")
                                            graph.add_node(bnd_name, bnd_name_a_g)

                                            sbnd1_name = gen_node("s_BND", sbnd1_name_a_g, bnd_name)
                                            sbnd2_name = gen_node("s_BND", sbnd2_name_a_g, bnd_name)

                                            action_graph.add_edge(a1_type+"-"+s1, sbnd1_name_a_g)
                                            action_graph.add_edge(a2_type+"-"+s2, sbnd2_name_a_g)
                                            graph.add_edge(sites_names[i][k], sbnd1_name)
                                            graph.add_edge(sites_names[j][l], sbnd2_name)
            elif a1 in RHS.keys():

                # a1 have been created

                for j in range(len(agents)):
                    if i != j :
                        a2 = agents[j]
                        if a2 in RHS.keys() and\
                           {a2, a1} not in bounds:

                            # a2 have been created or was already here

                            a2_type = RHS[a2]['name']

                            for k in range(len(a_sites[i])):
                                for l in range(k, len(a_sites[j])):
                                    s1 = a_sites[i][k][0]
                                    s2 = a_sites[j][l][0]

                                    if s1 in RHS[a1]['sites'].keys() and\
                                       s2 in RHS[a2]['sites'].keys():

                                        new_bind1 = RHS[a1]['sites'][s1]["binding"]
                                        new_bind2 = RHS[a2]['sites'][s2]["binding"]

                                        if new_bind1 == new_bind2 and (new_bind1 is not None):

                                                # BND between two syntesized nodes

                                                bnd_name_a_g = gen_name("BND", action_graph)
                                                action_graph.add_node(bnd_name_a_g, "BND")

                                                sbnd1_name_a_g = gen_node("s_BND", "s_BND", bnd_name_a_g, action_graph)
                                                sbnd2_name_a_g = gen_node("s_BND", "s_BND", bnd_name_a_g, action_graph)

                                                bnd_name = gen_name("BND")
                                                graph.add_node(bnd_name, bnd_name_a_g)

                                                sbnd1_name = gen_node("s_BND", sbnd1_name_a_g, bnd_name)
                                                sbnd2_name = gen_node("s_BND", sbnd2_name_a_g, bnd_name)

                                                action_graph.add_edge(a1_type+"-"+s1, sbnd1_name_a_g)
                                                action_graph.add_edge(a2_type+"-"+s2, sbnd2_name_a_g)
                                                graph.add_edge(sites_names[i][k], sbnd1_name)
                                                graph.add_edge(sites_names[j][l], sbnd2_name)

                                                bounds.append({a1,a2})


    @staticmethod
    def uncompile_rules(rules, action_graph):
        graph = TypedDiGraph()
        for i in range(len(rules)):
            KappaImporter.uncompile_rule(rules[i], graph, action_graph, str(i))
        return graph

    @staticmethod
    def uncompile(files, parser="kappa_to_graph.byte", out="out", del_out=True, list_=True):
        if out[-1] != "/":
            out += "/"

        if not os.path.exists(out):
            os.makedirs(out)

        if list_ :
            nug_list = []

            subprocess.check_call(("./%s %s -path %s" % (parser, " ".join(files), out)).split(" "))

            f_a = open(out+'agents.json', 'r+')
            f_r = open(out+'rules.json', 'r+')

            agents_ = json.loads(f_a.read())
            action_graph = KappaImporter.uncompile_agents(agents_)
            rules = json.loads(f_r.read())
            for rule in rules:
                nugget = TypedDiGraph()
                KappaImporter.uncompile_rule(rule, nugget, action_graph)

                nugget.metamodel_ = action_graph
                nugget.hom = TypedHomomorphism.canonic(nugget, action_graph)

                nug_list.append(nugget)

            if del_out:
                subprocess.check_call(("rm -rf %s %s %s" % (out+"rules.ka", out+"rules.json", out+"agents.json")).split(" "))
                try:
                    subprocess.check_call(("rmdir %s" % out).split(" "))
                except subprocess.CalledProcessError:
                    pass

            return nug_list
        else:
            subprocess.check_call(("./%s %s -path %s" % (parser, " ".join(files), out)).split(" "))

            f_a = open(out+'agents.json', 'r+')
            f_r = open(out+'rules.json', 'r+')

            agents = json.loads(f_a.read())
            rules = json.loads(f_r.read())

            action_graph = KappaImporter.uncompile_agents(agents)
            graph = KappaImporter.uncompile_rules(rules, action_graph)

            graph.metamodel_ = action_graph
            graph.hom = TypedHomomorphism.canonic(graph, action_graph)

            if del_out:
                subprocess.check_call(("rm -rf %s" % (out)).split(" "))

            return graph
