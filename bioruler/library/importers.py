"""."""
from regraph.library.data_structures import TypedDiGraph

from bioruler.library.metamodels import metamodel_AG
from bioruler.library.biopax_utils import BioPAXModel

# from bioruler.library.utils import (edge_from_ids,
#                                     generate_family_action_node,
#                                     generate_family_source,
#                                     generate_family_target)


def generate_family_action_node(family_id):
    return ("%s_family" % family_id, "FAM")


def generate_family_source(family_id):
    return ("%s_family_s" % family_id, "FAM_s")


def generate_family_target(family_id):
    return ("%s_family_t" % family_id, "FAM_t")


class BioPaxActionGraphImporter():
    """."""

    def __init__(self):
        """Initialize BioPAX model."""
        self.metamodel_ = metamodel_AG
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
                complexes_data[uri]["components"] = set()
                for el in complex.getComponent():
                    if not el.getModelInterface() == self.data_.complex_class_:
                        complexes_data[uri]["components"].add(el.getEntityReference().getUri())
                    else:
                        complexes_data[uri]["components"].add(el.getUri())
                # Collect residues and flags
                modifications = self.data_.get_modifications([uri])
                complexes_data[uri].update(modifications)
            else:
                if not ignore_families:
                    families[uri] = {}
                    families[uri]["members"] = set()
                    for m in complex.getMemberPhysicalEntity():
                        if not m.getModelInterface() == self.data_.complex_class_:
                            families[uri]["members"].add(m.getEntityReference().getUri())
                        else:
                            families[uri]["members"].add(m.getUri())

                    # Collect residues and flags
                    modifications = self.data_.get_modifications([uri])
                    families[uri].update(modifications)

        return (complexes_data, families)

    def collect_modifications(self, ignore_families=False):
        pass

    def generate_proteins(self, proteins, graph):
        nodes = []
        edges = []

        for protein_id, data in proteins.items():
            # Create nodes representing proteins
            protein_node = self.data_.protein_reference_to_node(protein_id)
            if protein_node is not None:
                nodes.append(protein_node)
                for flag in data["flags"]:
                        flag_node = self.data_.flag_to_node(flag, protein_id)
                        if flag_node is not None:
                            nodes.append(flag_node)

                            edges.append(
                                (flag_node[0], protein_id)
                            )
                for residue in data["residues"]:
                    residue_node = self.data_.residue_to_node(residue, protein_id)
                    if residue_node is not None:
                        nodes.append(residue_node)

                        edges.append(
                            (residue_node[0], protein_id)
                        )
                        flag_node = self.data_.flag_to_node(residue, residue_node[0])
                        if flag_node is not None:
                            nodes.append(flag_node)
                            edges.append(
                                (flag_node[0], residue_node[0])
                            )
                # Create nodes representing regions and region attributes
                for region_id, region_data in data["regions"].items():
                    region_node = self.data_.region_to_node(region_id)
                    if region_node is not None:
                        nodes.append(region_node)
                        edges.append(
                            (region_id, protein_id)
                        )
                        for flag in region_data["flags"]:
                            flag_node = self.data_.flag_to_node(flag, region_id)
                            if flag_node is not None:
                                nodes.append(flag_node)
                                edges.append(
                                    (flag_node[0], region_id)
                                )
                        for residue in region_data["residues"]:
                            residue_node = self.data_.residue_to_node(residue, region_id)
                            if residue_node is not None:
                                nodes.append(residue_node)
                                edges.append(
                                    (residue_node[0], region_id)
                                )
                                edges.append(
                                    (residue_node[0], protein_id)
                                )

                                flag_node = self.data_.flag_to_node(residue, residue_node[0])
                                if flag_node is not None:
                                    nodes.append(flag_node)
                                    edges.append(
                                        (flag_node[0], residue_node[0])
                                    )
        graph.add_nodes_from(nodes)
        graph.add_edges_from(edges)

    def generate_families(self, families, graph):
        nodes = []
        edges = []
        for family_id, data in families.items():
            # Create family node
            family_node = self.data_.family_reference_to_node(family_id)
            if family_node is not None:
                nodes.append(family_node)
                # Create flag nodes
                if "flag" in data.keys():
                    for flag in data["flags"]:
                            flag_node = self.data_.flag_to_node(flag, family_id)
                            if flag_node is not None:
                                nodes.append(flag_node)
                                edges.append(
                                    (flag_node[0], family_node[0])
                                )
                # Create residue nodes
                if "residue" in data.keys():
                    for residue in data["residues"]:
                        residue_node = self.residue_to_node(residue, family_id)
                        if residue_node is not None:
                            nodes.append(residue_node)
                            edges.append(
                                (residue_node[0], family_node[0])
                            )
                            flag_node = self.data_.flag_to_node(
                                residue, residue_node[0])
                            if flag_node is not None:
                                nodes.append(flag_node)
                                edges.append(
                                    (flag_node[0], residue_node[0])
                                )
                # Create FAM action
                fam = generate_family_action_node(family_id)
                fam_s = generate_family_source(family_id)
                fam_t = generate_family_target(family_id)

                nodes.append(fam)
                nodes.append(fam_s)
                nodes.append(fam_t)

                edges.append((fam_s[0], fam[0]))
                edges.append((fam_t[0], fam[0]))

                edges.append((fam_t[0], family_node[0]))

                for member in data["members"]:
                    if member in graph.nodes():
                        edges.append((member, fam_s[0]))

        graph.add_nodes_from(nodes)
        graph.add_edges_from(edges)

    def generate_small_molecules(self, small_molecules, graph):
        nodes = []
        for molecule_id in small_molecules:
            molecule_node = self.data_.small_molecule_to_node(molecule_id)
            if molecule_node is not None:
                nodes.append(molecule_node)
        graph.add_nodes_from(nodes)

    def generate_complexes(self, complexes, graph):
        nodes = []
        edges = []

        for complex_id, data in complexes.items():
            complex_node = self.data_.complex_to_node(complex_id)
            nodes.append(complex_node)

            # Create flag nodes
            if "flag" in data.keys():
                for flag in data["flags"]:
                        flag_node = self.data_.flag_to_node(flag, complex_id)
                        if flag_node is not None:
                            nodes.append(flag_node)
                            edges.append(
                                (flag_node[0], complex_node[0])
                            )
            # Create residue nodes
            if "residue" in data.keys():
                for residue in data["residues"]:
                    residue_node = self.residue_to_node(residue, complex_id)
                    if residue_node is not None:
                        nodes.append(residue_node)
                        edges.append(
                            (residue_node[0], complex_node[0])
                        )
                        flag_node = self.data_.flag_to_node(residue, residue_node['id'])
                        if flag_node is not None:
                            nodes.append(flag_node)
                            edges.append(
                                (flag_node[0], residue_node[0])
                            )

            for component in data["components"]:
                if component in graph.nodes():
                    edges.append(
                        (component, complex_node[0])
                    )

        graph.add_nodes_from(nodes)
        graph.add_edges_from(edges)

    # def complex(self, graph):
    #     print("GENERATING COMPLEXES...")
    #     nodes = []

    #     if self.data_ is not None:
    #         complexes = self.data_.model_.getObjects(self.data_.complex_class_)
    #     else:
    #         raise ValueError("BioPAX model is not loaded!")

    #     for complex in list(complexes):
    #         uri = complex.getUri()
    #         if not self.data_.is_complex_family(uri):
    #             complex_node = self.data_.complex_to_node(uri)
    #             if complex_node is not None:
    #                 nodes.append(complex_node)
    #     graph.add_nodes_from(nodes)

    # def protein(self, graph):
    #     print("GENERATING PROTEINS...")
    #     nodes = []

    #     if self.data_ is not None:
    #         proteins = self.data_.model_.getObjects(self.data_.protein_reference_class_)
    #     else:
    #         raise ValueError("BioPAX model is not loaded!")

    #     # Iterate through protein references
    #     for protein in list(proteins):
    #         uri = protein.getUri()
    #         if not self.data_.is_protein_family(uri):
    #             protein_node = self.data_.protein_reference_to_node(uri)
    #             if protein_node is not None:
    #                 nodes.append(protein_node)
    #     graph.add_nodes_from(nodes)

    # def region(self, graph):
    #     print("GENERATING REGIONS...")
    #     nodes = []

    #     if self.data_ is not None:
    #         proteins = self.data_.model_.getObjects(self.data_.protein_class_)
    #     else:
    #         raise ValueError("BioPAX model is not loaded!")

    #     region_refrences = set()
    #     for protein in proteins:
    #         if self.data_.is_fragment(protein.getUri()):
    #             features = protein.getFeature()
    #             for f in features:
    #                 if f.getModelInterface() == self.data_.fragment_feature_class_:
    #                     reference_proteins = set(
    #                         [el.getEntityReference().getUri() for el in f.getFeatureOf()]
    #                     )
    #                     if len(reference_proteins) > 1:
    #                         raise ValueError("Region references to more than one protein!")
    #                     region_refrences.add(f.getUri())

    #     for region in region_refrences:
    #         region_node = self.data_.region_to_node(region)
    #         if region_node is not None:
    #             nodes.append(region_node)

    #     graph.add_nodes_from(nodes)

    # def residue(self, graph):
    #     print("GENRATING RESIDUE...")
    #     nodes = []
    #     # in our metamodel residue can belong to the
    #     # - protein
    #     # - complex

    #     if self.data_ is not None:
    #         proteins = self.data_.model_.getObjects(self.data_.protein_class_)
    #         complexes = self.data_.model_.getObjects(self.data_.complex_class_)
    #     else:
    #         raise ValueError("BioPAX model is not loaded!")

    #     residues = set()
    #     for protein in proteins:
    #         for f in protein.feature:
    #             if f.getModelInterface() == self.data_.modification_feature_class_:
    #                 if f.getFeatureLocation() is not None:
    #                     residues.add(f.getUri())
    #     for complex in complexes:
    #         for f in complex.feature:
    #             if f.getModelInterface() == self.data_.modification_feature_class_:
    #                 if f.getFeatureLocation() is not None:
    #                     residues.add(f.getUri())
    #     for residue in residues:
    #         residue_node = self.data_.residue_to_node(residue)
    #         if residue_node is not None:
    #             nodes.append(residue_node)
    #     graph.add_nodes_from(nodes)

    # def state(self, graph):
    #     print("GENERATING SMALL MOLECULES...")
    #     nodes = []

    #     if self.data_ is not None:
    #         proteins = self.data_.model_.getObjects(self.data_.protein_class_)
    #         complexes = self.data_.model_.getObjects(self.data_.complex_class_)
    #         small_molecules = self.data_.model_.getObjects(self.data_.small_molecule_class_)
    #     else:
    #         raise ValueError("BioPAX model is not loaded!")

    #     flags = set()
    #     for protein in proteins:
    #         for f in protein.feature:
    #             if f.getModelInterface() == self.data_.modification_feature_class_:
    #                 if f.getFeatureLocation() is None:
    #                     flags.add(f.getUri())
    #     for complex in complexes:
    #         for f in complex.feature:
    #             if f.getModelInterface() == self.data_.modification_feature_class_:
    #                 if f.getFeatureLocation() is None:
    #                     flags.add(f.getUri())
    #     for molecule in small_molecules:
    #         for f in molecule.feature:
    #             if f.getModelInterface() == self.data_.modification_feature_class_:
    #                 if f.getFeatureLocation() is None:
    #                     flags.add(f.getUri())
    #     for flag in flags:
    #         flag_node = self.data_.flag_to_node(flag)
    #         if flag_node is not None:
    #             nodes.append(flag_node)
    #     graph.add_nodes_from(nodes)

    # def small_molecule(self, graph):
    #     print("GENERATING flags...")
    #     # protein flag
    #     # complex flag
    #     # small molecule flag
    #     nodes = []

    #     if self.data_ is not None:
    #         small_molecules = self.data_.model_.getObjects(self.data_.small_molecule_reference_class_)
    #     else:
    #         raise ValueError("BioPAX model is not loaded!")

    #     # Iterate through protein references
    #     for molecule in list(small_molecules):
    #         uri = molecule.getUri()
    #         if not self.data_.is_protein_family(uri):
    #             molecule_node = self.data_.small_molecule_to_node(uri)
    #             if molecule_node is not None:
    #                 nodes.append(molecule_node)
    #     graph.add_nodes_from(nodes)

    # def family(self, graph):
    #     print("GENERATING Family...")
    #     nodes = []

    #     if self.data_ is not None:
    #         proteins = self.data_.model_.getObjects(self.data_.protein_reference_class_)
    #         complexes = self.data_.model_.getObjects(self.data_.complex_class_)
    #         small_molecules = self.data_.model_.getObjects(self.data_.small_molecule_reference_class_)
    #     else:
    #         raise ValueError("BioPAX model is not loaded!")
    #     for protein in proteins:
    #         if self.data_.is_protein_family(protein.getUri()):
    #             family_node = self.data_.family_reference_to_node(protein.getUri())
    #             if family_node is not None:
    #                 nodes.append(family_node)
    #     for complex in complexes:
    #         if self.data_.is_complex_family(complex.getUri()):
    #             family_node = self.data_.family_reference_to_node(complex.getUri())
    #             if family_node is not None:
    #                 nodes.append(family_node)
    #     for molecule in small_molecules:
    #         if self.data_.is_protein_family(molecule.getUri()):
    #             family_node = self.data_.family_reference_to_node(molecule.getUri())
    #             if family_node is not None:
    #                 nodes.append(family_node)

    #     graph.add_nodes_from(nodes)

    # def BND(self, graph):
    #     print("GENRATING BND")

    # def BND_s(self, graph):
    #     print("GENRATING BND_s")

    # def BRK(self, graph):
    #     print("GENRATING BRK")

    # def BRK_t(self, graph):
    #     print("GENRATING BRK_t")

    # def MOD(self, graph):
    #     print("GENRATING MOD")

    # def MOD_s(self, graph):
    #     print("GENRATING MOD_s")

    # def MOD_t(self, graph):
    #     print("GENRATING MOD_t")

    # def FAM(self, graph):
    #     print("GENRATING FAM")

    # def FAM_s(self, graph):
    #     print("GENRATING FAM_s")

    # def FAM_t(self, graph):
    #     print("GENRATING FAM_t")

    def collect_agents(self, ignore_families=False):
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

    def collect_actions(self, ignore_families=False):
        modifications = self.collect_modifications(ignore_families)
        actions = {
            "MOD": modifications
        }
        return actions

    def import_model(self, filename, ignore_families=False):
        self.data_.load(filename)
        graph = TypedDiGraph(self.metamodel_)
        # for node_type in self.metamodel_.nodes():
        #     if self.metamodel_.node[node_type].type_ != "action":
        #         getattr(self, node_type)(graph)
        #         # print(graph.nodes())
        #         print("# nodes: ", len(graph.nodes()))
        #         print(self.metamodel_.neighbors(node_type))

        # return graph

        agents = self.collect_agents(ignore_families)
        actions = self.collect_actions(ignore_families)

        self.generate_proteins(agents["proteins"], graph)
        self.generate_families(agents["protein_families"], graph)
        self.generate_small_molecules(agents["small_molecules"], graph)
        self.generate_families(agents["small_molecule_families"], graph)
        self.generate_complexes(agents["complexes"], graph)
        self.generate_families(agents["complex_families"], graph)

        return graph
