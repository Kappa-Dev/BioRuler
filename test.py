"""."""
from bioruler.library.importers import BioPAXImporter


if __name__ == '__main__':
    a = BioPAXImporter()
    graph = a.import_model("NCI-Nature-Curated-final-1.bp3.owl")
    print("Nodes: ", len(graph.nodes()))
    print("Edges: ", len(graph.edges()))
