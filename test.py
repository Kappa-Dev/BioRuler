"""."""
from bioruler.library.importers import BioPaxActionGraphImporter


if __name__ == '__main__':
    a = BioPaxActionGraphImporter()
    graph = a.import_model("NCI-Nature-Curated-final-1.bp3.owl")
    print("Nodes: ", len(graph.nodes()))
    print("Edges: ", len(graph.edges()))
    graph.export("Graph.json")
