"""."""
from bioruler.library.importers import BioPAXImporter

from regraph.library.utils import plot_graph


if __name__ == '__main__':
    a = BioPAXImporter()
    # Test the importer on PID database (https://pid.nci.nih.gov/)
    # the result of an import is the action graph and the collection
    # of nuggets. Mapping between the nodes of the action graph and the
    # nuggets is implicit (through the same ids of the nodes)
    action_graph, nuggets = a.import_model("data/NCI-Nature-Curated-final-1.bp3.owl")
    print("Action graph nodes: ", len(action_graph.nodes()))
    print("Action graph edges: ", len(action_graph.edges()))
    print("Number of nuggets: ", len(nuggets))

    # print a couple of nuggets
    plot_graph(nuggets[0], filename="nugget1.png")
    print("NUGGET 1\n")
    for n in nuggets[0]:
        print(nuggets[0].node[n])
    print()

    plot_graph(nuggets[10], filename="nugget2.png")
    print("NUGGET 2\n")
    for n in nuggets[10]:
        print(nuggets[10].node[n])
    print()

    plot_graph(nuggets[100], filename="nugget3.png")
    print("NUGGET 3\n")
    for n in nuggets[100]:
        print(nuggets[100].node[n])
    print()

    # We can export graph to json
    action_graph.export("data/Graph.json")
