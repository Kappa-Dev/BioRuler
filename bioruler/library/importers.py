"""."""
import json

from regraph.library.data_structures import TypedDiGraph

from bioruler.library.metamodels import simple_metamodel_AG
from bioruler.library.biopax_utils import (BioPAXModel, edge_from_ids,
                                           generate_family_action_node)


class BioPAXImporter():
    """."""

    def __init__(self):
        """Initialize BioPAX model."""
        self.data_ = BioPAXModel()

    def collect_proteins(self, ignore_families=False):
        """."""
        if self.data_ is not None:
            proteins = self.data_.model_.getObjects(self.data_.protein_reference_class_)
        else:
            raise ValueError("BioPAX model is not loaded!")
        protein_data = {}
        families = {}

        # Iterate through protein references
        for protein in list(proteins):
            uri = protein.getUri()

            # Individual proteins
            if not self.data_.is_protein_family(protein.getUri()):
                protein_data[uri] = {}

                physical_proteins = protein.getEntityReferenceOf()
                distinct_regions = {}
                for el in physical_proteins:
                    # Collect regions
                    if self.data_.is_fragment(el.getUri()):
                        features = el.getFeature()
                        for f in features:
                            if f.getModelInterface() == self.data_.fragment_feature_class_:
                                if f.getUri() in distinct_regions.keys():
                                    distinct_regions[f.getUri()].add(
                                        el.getUri())
                                else:
                                    distinct_regions[f.getUri()] = set()
                regions = {}
                # Get PTM for each region
                for region_id, physical_entities_id in distinct_regions.items():
                    regions[region_id] = self.data_.get_modifications(
                        physical_entities_id)
                protein_data[uri]["regions"] = regions

                # Collect residues and flags
                physical_proteins_id =\
                    [el.getUri() for el in physical_proteins]
                residues = set()
                modifications = self.data_.get_modifications(physical_proteins_id)
                # If residue is in the location of region, we add it to the region
                for region_id in distinct_regions.keys():
                    for residue_id in modifications["residues"]:
                        if self.data_.residue_in_region(residue_id, region_id):
                            regions[region_id]["residues"].add(residue_id)
                        else:
                            residues.add(residue_id)
                protein_data[uri]["residues"] = residues
                protein_data[uri]["flags"] = modifications["flags"]
            # Families
            else:
                if not ignore_families:
                    if len(protein.getMemberEntityReference()) > 0:
                        families[uri] = {}
                        families[uri]["members"] =\
                            set([m.getUri() for m in protein.getMemberEntityReference()])

                        physical_proteins = protein.getEntityReferenceOf()

                        # Collect residues and flags
                        physical_proteins_id = [el.getUri() for el in physical_proteins]
                        modifications = self.data_.get_modifications(physical_proteins_id)
                        families[uri].update(modifications)
                    else:
                        families[uri] = {}
                        physical_entities = protein.getEntityReferenceOf()
                        for entity in physical_entities:
                            if len(entity.getMemberPhysicalEntity()) > 0:
                                if "members" in families[uri].keys():
                                    families[uri]["members"].update(
                                        [m.getEntityReference().getUri() for m in entity.getMemberPhysicalEntity()])
                                else:
                                    families[uri]["members"] =\
                                        set([m.getEntityReference().getUri() for m in entity.getMemberPhysicalEntity()])
                        physical_entities_ids = [m.getUri() for m in physical_entities]
                        modifications = self.data_.get_modifications(physical_entities_ids)
                        families[uri].update(modifications)
        return (protein_data, families)

    def collect_small_molecules(self, ignore_families=False):
        small_molecules = self.data_.model_.getObjects(
            self.data_.small_molecule_reference_class_)
        small_molecules_data = set()
        families = {}
        for molecule in small_molecules:
            uri = molecule.getUri()
            if not self.data_.is_protein_family(molecule.getUri()):
                small_molecules_data.add(uri)
            else:
                if not ignore_families:
                    if len(molecule.getMemberEntityReference()) > 0:
                        families[uri] = {}
                        families[uri]["members"] =\
                            set([m.getUri() for m in molecule.getMemberEntityReference()])
                    else:
                        families[uri] = {}
                        physical_entities = molecule.getEntityReferenceOf()
                        for entity in physical_entities:
                            if len(entity.getMemberPhysicalEntity()) > 0:
                                if "members" in families[uri].keys():
                                    families[uri]["members"].update(
                                        [m.getEntityReference().getUri() for m in entity.getMemberPhysicalEntity()])
                                else:
                                    families[uri]["members"] =\
                                        set([m.getEntityReference().getUri() for m in entity.getMemberPhysicalEntity()])
        return (small_molecules_data, families)

    def collect_complexes(self, ignore_families=False):
        complexes = self.data_.model_.getObjects(
            self.data_.complex_class_)

        complexes_data = {}
        families = {}

        for complex in complexes:
            uri = complex.getUri()
            if not self.data_.is_complex_family(uri):
                complexes_data[uri] = {}
                complexes_data[uri]["components"] =\
                    set([el.getUri() for el in complex.getComponent()])
                # Collect residues and flags
                modifications = self.data_.get_modifications([uri])
                complexes_data[uri].update(modifications)
            else:
                if not ignore_families:
                    families[uri] = {}
                    families[uri]["members"] =\
                        set([m.getUri() for m in complex.getMemberPhysicalEntity()])

                    # Collect residues and flags
                    modifications = self.data_.get_modifications([uri])
                    families[uri].update(modifications)

        return (complexes_data, families)

    def collect_agents_simple_mm_AG(self, ignore_families=False):
        (proteins, protein_families) = self.collect_proteins(ignore_families)
        (small_molecules, small_molecule_families) =\
            self.collect_small_molecules(ignore_families)
        (complexes, complex_families) =\
            self.collect_complexes(ignore_families)
        agents = {
            "proteins": proteins,
            "protein_families": protein_families,
            "small_molecules": small_molecules,
            "small_molecule_families": small_molecule_families,
            "complexes": complexes,
            "complex_families": complex_families
        }
        return agents

    def convert_collection_protein(self, proteins):
        agent_nodes = []
        structural_edges = []

        for protein_id, data in proteins.items():
            # Create nodes representing proteins
            protein_node = self.data_.protein_reference_to_node(protein_id)
            if protein_node is not None:
                agent_nodes.append(protein_node)
                for flag in data["flags"]:
                        flag_node = self.data_.flag_to_node(protein_id, flag)
                        if flag_node is not None:
                            agent_nodes.append(flag_node)

                            structural_edges.append(
                                edge_from_ids(
                                    flag_node["id"],
                                    protein_node["id"]
                                )
                            )
                for residue in data["residues"]:
                    residue_node = self.data_.residue_to_node(protein_id, residue)
                    if residue_node is not None:
                        agent_nodes.append(residue_node)

                        structural_edges.append(
                            edge_from_ids(
                                residue_node["id"],
                                protein_node["id"]
                            )
                        )
                        flag_node = self.data_.flag_to_node(residue_node['id'], residue)
                        if flag_node is not None:
                            agent_nodes.append(flag_node)
                            structural_edges.append(
                                edge_from_ids(
                                    flag_node["id"],
                                    residue_node["id"]
                                )
                            )
                # Create nodes representing regions and region attributes
                for region_id, region_data in data["regions"].items():
                    region_node = self.data_.region_to_node(region_id)
                    if region_node is not None:
                        agent_nodes.append(region_node)
                        structural_edges.append(
                            edge_from_ids(
                                region_node["id"],
                                protein_node["id"]
                            )
                        )
                        for flag in region_data["flags"]:
                            flag_node = self.data_.flag_to_node(region_id, flag)
                            if flag_node is not None:
                                agent_nodes.append(flag_node)
                                structural_edges.append(
                                    edge_from_ids(
                                        flag_node["id"],
                                        region_node["id"]
                                    )
                                )
                        for residue in region_data["residues"]:
                            residue_node = self.data_.residue_to_node(region_id, residue)
                            if residue_node is not None:
                                agent_nodes.append(residue_node)
                                structural_edges.append(
                                    edge_from_ids(
                                        residue_node["id"],
                                        region_node["id"]
                                    )
                                )
                                structural_edges.append(
                                    edge_from_ids(
                                        residue_node["id"],
                                        protein_node["id"]
                                    )
                                )

                                flag_node = self.data_.flag_to_node(residue_node["id"], flag)
                                if flag_node is not None:
                                    agent_nodes.append(flag_node)
                                    structural_edges.append(
                                        edge_from_ids(
                                            flag_node["id"],
                                            residue_node["id"]
                                        )
                                    )
        return (agent_nodes, structural_edges)

    def convert_collection_family(self, families):
        agent_nodes = []
        structural_edges = []
        for family_id, data in families.items():
            # Create family node
            family_node = self.data_.family_reference_to_node(family_id)
            if family_node is not None:
                agent_nodes.append(family_node)
                # Create flag nodes
                if "flag" in data.keys():
                    for flag in data["flags"]:
                            flag_node = self.data_.flag_to_node(family_id, flag)
                            if flag_node is not None:
                                agent_nodes.append(flag_node)
                                structural_edges.append(
                                    edge_from_ids(
                                        flag_node["id"],
                                        family_node["id"]
                                    )
                                )
                # Create residue nodes
                if "residue" in data.keys():
                    for residue in data["residues"]:
                        residue_node = self.residue_to_node(family_id, residue)
                        if residue_node is not None:
                            agent_nodes.append(residue_node)
                            structural_edges.append(
                                edge_from_ids(
                                    residue_node["id"],
                                    family_node["id"]
                                )
                            )
                            flag_node = self.data_.flag_to_node(residue_node['id'], residue)
                            if flag_node is not None:
                                agent_nodes.append(flag_node)
                                structural_edges.append(
                                    edge_from_ids(
                                        flag_node["id"],
                                        residue_node["id"]
                                    )
                                )
                # Create FAM action
                action_node = generate_family_action_node(family_id)
                agent_nodes.append(action_node)
                structural_edges.append(
                    edge_from_ids(
                        action_node["id"],
                        family_node["id"]
                    )
                )
                for member in data["members"]:
                    structural_edges.append(
                        edge_from_ids(
                            member,
                            action_node["id"]
                        )
                    )
        return (agent_nodes, structural_edges)

    def convert_collection_small_molecule(self, small_molecules):
        agent_nodes = []
        structural_edges = []
        for molecule_id in small_molecules:
            molecule_node = self.data_.small_molecule_to_node(molecule_id)
            if molecule_node is not None:
                agent_nodes.append(molecule_node)
        return (agent_nodes, structural_edges)

    def convert_collection_complex(self, small_molecules):
        pass

    def import_model(self, filename, metamodel=None, ignore_families=False):
        if metamodel is None:
            metamodel = simple_metamodel_AG
        self.data_.load(filename)
        agents = self.collect_agents_simple_mm_AG(ignore_families)
        # print(agents["protein_families"])
        graph = {"nodes": [], "edges": []}
        nodes, edges = self.convert_collection_protein(agents["proteins"])
        graph["nodes"] += nodes
        graph["edges"] += edges
        nodes, edges = self.convert_collection_family(agents["protein_families"])
        graph["nodes"] += nodes
        graph["edges"] += edges
        nodes, edges = self.convert_collection_small_molecule(agents["small_molecules"])
        graph["nodes"] += nodes
        graph["edges"] += edges
        nodes, edges = self.convert_collection_family(agents["small_molecule_families"])
        graph["nodes"] += nodes
        graph["edges"] += edges

        with open("BioPAX_gaph.json", "w") as f:
            json.dump(graph, f, sort_keys=True, indent=4)

        res_graph = TypedDiGraph(simple_metamodel_AG)
        res_graph.load("BioPAX_gaph.json")
        return res_graph
        # for node in metamodel:
            # print(node)
            # print('%s_to_graph' % node)
            # getattr(self, '%s_to_graph' % node)()
        # actions = self.collect_actions_simple_mm_AG()
        # print(agents)

