"""The importers from the biological data fromats to the action graph and nuggets."""
from regraph.library.data_structures import TypedDiGraph

from bioruler.library.metamodels import metamodel_AG
from bioruler.library.biopax_utils import BioPAXModel


def generate_family_action_node(family_id):
    """Generate FAM action node tuple."""
    return ("%s_family" % family_id, "FAM")


def generate_family_source(family_id):
    """Generate source FAM action node tuple."""
    return ("%s_family_s" % family_id, "FAM_s")


def generate_family_target(family_id):
    """Generate target FAM action node tuple."""
    return ("%s_family_t" % family_id, "FAM_t")


def generate_modification_source(reaction_id, target):
    """Generate source MOD action node tuple."""
    return ("%s_of_%s_s" % (reaction_id, target), "MOD_s")


def generate_modification_target(reaction_id, target):
    """Generate target MOD action node tuple."""
    return ("%s_of_%s_t" % (reaction_id, target), "MOD_t")


class BioPAXImporter():
    """."""

    def __init__(self):
        """Initialize BioPAX model."""
        self.metamodel_ = metamodel_AG
        self.data_ = BioPAXModel()

    def collect_proteins(self, ignore_families=False):
        """Collect the ids of proteins, regions, flags and residues."""
        print("Collecting proteins and families...")
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
                            if f.getModelInterface() ==\
                               self.data_.fragment_feature_class_:
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
                # Collect residues and flags
                physical_proteins_id =\
                    [el.getUri() for el in physical_proteins]
                modifications = self.data_.get_modifications(
                    physical_proteins_id)
                # If residue is in the location of the region, we add it to the region

                for region_id in distinct_regions.keys():
                    for residue_id in modifications["residues"]:
                        if self.data_.residue_in_region(residue_id, region_id):
                            regions[region_id]["residues"].add(residue_id)

                protein_data[uri]["residues"] = modifications["residues"]
                protein_data[uri]["flags"] = modifications["flags"]
                protein_data[uri]["regions"] = regions

            # Families
            else:
                if not ignore_families:
                    if len(protein.getMemberEntityReference()) > 0:
                        families[uri] = {}
                        families[uri]["members"] =\
                            set([m.getUri() for m in protein.getMemberEntityReference()])

                        physical_proteins = protein.getEntityReferenceOf()

                        # Collect residues and flags
                        physical_proteins_id =\
                            [el.getUri() for el in physical_proteins]
                        modifications = self.data_.get_modifications(
                            physical_proteins_id)
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
        """Collect ids of small molecules."""
        print("Collecting small molecules and families...")
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
        """Collect ids of complexes, flags, residues and components."""
        print("Collecting complexes and families...")
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
                        complexes_data[uri]["components"].add(
                            el.getEntityReference().getUri())
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
                            families[uri]["members"].add(
                                m.getEntityReference().getUri())
                        else:
                            families[uri]["members"].add(m.getUri())

                    # Collect residues and flags
                    modifications = self.data_.get_modifications([uri])
                    families[uri].update(modifications)

        return (complexes_data, families)

    def collect_modifications(self, ignore_families=False):
        """Detect and collect modifications and their participants."""
        print("Collecting modifications...")
        catalysis = self.data_.model_.getObjects(
            self.data_.catalysis_class_)
        modification_data = {}
        for reaction in catalysis:
            uri = reaction.getUri()
            # Get reaction controlled
            for controlled_reaction in reaction.getControlled():
                # We are intereseted only in Biochemical Reactions
                if controlled_reaction.getModelInterface() ==\
                   self.data_.biochemical_reaction_class_:
                    LHS = controlled_reaction.getLeft()
                    RHS = controlled_reaction.getRight()

                    # Extract trivial modifications

                    # Check if it is a modification of the single
                    # entity's state
                    if len(LHS) == 1 and len(RHS) == 1:
                        initial_entity = list(LHS)[0]
                        resulting_entity = list(RHS)[0]

                        elements_match = False
                        if initial_entity.getModelInterface() ==\
                           resulting_entity.getModelInterface():
                            if initial_entity.getModelInterface() !=\
                               self.data_.complex_class_:
                                if initial_entity.getEntityReference().getUri() ==\
                                   resulting_entity.getEntityReference().getUri():
                                    elements_match = True
                            else:
                                if initial_entity.getUri() == resulting_entity.getUri():
                                    elements_match = True
                        # If both sides have the same entity
                        if elements_match:
                            modification_data[uri] =\
                                {"sources": {}, "targets": []}
                            for controller in reaction.getController():
                                if controller.getModelInterface() !=\
                                   self.data_.complex_class_:
                                    entity = controller.getEntityReference().getUri()
                                else:
                                    entity = controller.getUri()
                                modification_data[uri]["sources"][entity] = []
                                for f in controller.getFeature():
                                    modification_data[uri]["sources"][entity].append(
                                        f.getUri()
                                    )

                            for f in resulting_entity.getFeature():
                                if resulting_entity.getModelInterface() !=\
                                   self.data_.complex_class_:
                                    if f not in initial_entity.getFeature():
                                        modification_data[uri]["targets"].append(
                                            (resulting_entity.getEntityReference().getUri(),
                                             f.getUri(), 1))
                                else:
                                    if f not in initial_entity.getFeature():
                                        modification_data[uri]["targets"].append(
                                            (resulting_entity.getUri(),
                                             f.getUri(), 1))
                            for f in initial_entity.getFeature():
                                if initial_entity.getModelInterface() !=\
                                   self.data_.complex_class_:
                                    if f not in resulting_entity.getFeature():
                                        modification_data[uri]["targets"].append(
                                            (resulting_entity.getEntityReference().getUri(),
                                             f.getUri(), 0))

                    # Extract phenomenological modifications (Voodoo starts here)
                    lhs_complexes = []
                    rhs_complexes = []
                    for el in LHS:
                        if el.getModelInterface() ==\
                           self.data_.complex_class_:
                            lhs_complexes.append(el)
                    for el in RHS:
                        if el.getModelInterface() ==\
                           self.data_.complex_class_:
                            rhs_complexes.append(el)

                    if len(lhs_complexes) == 1 and len(rhs_complexes) == 1:
                        if lhs_complexes[0].getUri() != rhs_complexes[0].getUri():
                            lhs_components = set([c.getUri() for c in lhs_complexes[0].getComponent()])
                            rhs_components = set([c.getUri() for c in rhs_complexes[0].getComponent()])
                            if lhs_components != rhs_components:
                                modification_data[uri] =\
                                    {"sources": {}, "targets": []}

                                for controller in reaction.getController():
                                    print(controller.getName())
                                    if controller.getModelInterface() !=\
                                       self.data_.complex_class_:
                                        entity = controller.getEntityReference().getUri()
                                    else:
                                        entity = controller.getUri()
                                    modification_data[uri]["sources"][entity] = []
                                    for f in controller.getFeature():
                                        print("\t\t", f)
                                        modification_data[uri]["sources"][entity].append(
                                            f.getUri()
                                        )

                                print("LHS ", lhs_complexes[0].getName())
                                for f in lhs_complexes[0].getFeature():
                                    print("\t", f)
                                for component in lhs_components:
                                    entity = self.data_.model_.getByID(component)
                                    print("\t", entity.getName())
                                    for f in entity.getFeature():
                                        print("\t\t", f)

                                print("RHS ", rhs_complexes[0].getName())
                                for f in rhs_complexes[0].getFeature():
                                    print("\t", f)
                                for component in rhs_components:
                                    entity = self.data_.model_.getByID(component)
                                    print("\t", entity.getName())
                                    for f in entity.getFeature():
                                        print("\t\t", f)

                                intersection = lhs_components.intersection(rhs_components)
                                l_difference = lhs_components.difference(rhs_components)
                                r_difference = rhs_components.difference(lhs_components)

                                l_ref = set()
                                left_ref_entity = []
                                for el in l_difference:
                                    entity = self.data_.model_.getByID(el)
                                    if entity.getModelInterface() != self.data_.complex_class_:
                                        l_ref.add(entity.getEntityReference().getUri())
                                        left_ref_entity.append((entity.getEntityReference().getUri(), entity))
                                    else:
                                        l_ref.add(entity.getUri())
                                        left_ref_entity.append((entity.getUri(), entity))

                                r_ref = set()
                                right_ref_entity = []
                                for el in r_difference:
                                    entity = self.data_.model_.getByID(el)
                                    if entity.getModelInterface() != self.data_.complex_class_:
                                        r_ref.add(entity.getEntityReference().getUri())
                                        right_ref_entity.append((entity.getEntityReference().getUri(), entity))
                                    else:
                                        r_ref.add(entity.getUri())
                                        right_ref_entity.append((entity.getUri(), entity))
                                if l_ref == r_ref:
                                    for ref, entity in right_ref_entity:
                                        print("\t", entity.getName())
                                        for left_ref, left_entity in left_ref_entity:
                                            if ref == left_ref:
                                                for f in entity.getFeature():
                                                    if f not in left_entity.getFeature():
                                                        print("\t\t\t", f, "1 - yes")
                                                        modification_data[uri]["targets"].append(
                                                            (ref, f.getUri(), 1)
                                                        )
                                                    else:
                                                        pass
                                                        print("\t\t\t", f, "1 - no")
                                                for f in left_entity.getFeature():
                                                    if f not in entity.getFeature():
                                                        print("\t\t\t", f, "0 - yes")
                                                        modification_data[uri]["targets"].append(
                                                            (ref, f.getUri(), 0)
                                                        )
                                                    else:
                                                        pass
                                                        print("\t\t\t", f, "0 - no")
                                else:
                                    if len(l_ref) == 1 and len(r_ref) == 1:
                                        left_entity = self.data_.model_.getByID(list(l_ref)[0])
                                        right_entity = self.data_.model_.getByID(list(r_ref)[0]) 
                                        if left_entity.getModelInterface() ==\
                                           self.data_.small_molecule_reference_class_ and\
                                           right_entity.getModelInterface() ==\
                                           self.data_.small_molecule_reference_class_:
                                            print("SMALL MOLECULE MODIFICATION: ")
                                            for el in intersection:
                                                entity = self.data_.model_.getByID(el)
                                                for f in rhs_complexes[0].getFeature():
                                                    if f not in lhs_complexes[0].getFeature():
                                                        print("\t\t\t", f, "1 - yes")
                                                        modification_data[uri]["targets"].append(
                                                            (entity.getEntityReference().getUri(), f.getUri(), 1)
                                                        )
                                                    else:
                                                        pass
                                                        print("\t\t\t", f, "1 - no")
                                                for f in lhs_complexes[0].getFeature():
                                                    if f not in rhs_complexes[0].getFeature():
                                                        print("\t\t\t", f, "0 - yes")
                                                        modification_data[uri]["targets"].append(
                                                            (entity.getEntityReference().getUri(), f.getUri(), 0)
                                                        )
                                                    else:
                                                        pass
                                                        print("\t\t\t", f, "0 - no")
        return modification_data

    # def collect_bindings(self, ignore_families=False):
    #     catalysis = self.data_.model_.getObjects(
    #         self.data_.catalysis_class_)
    #     modification_data = {}
    #     for reaction in catalysis:
    #         uri = reaction.getUri()
    #         # Get reaction controlled
    #         for controlled_reaction in reaction.getControlled():
    #             # We are intereseted in Complex Assemblies
    #             if controlled_reaction.getModelInterface() ==\
    #                self.data_.biochemical_reaction_class_:
    #                 LHS = controlled_reaction.getLeft()
    #                 RHS = controlled_reaction.getRight()

    def generate_proteins(self, proteins, graph):
        """Generate the nodes in the action graph for proteins."""
        print("Generating nodes for proteins...")
        nodes = []
        edges = []

        for protein_id, data in proteins.items():
            # Create nodes representing proteins
            protein_node = self.data_.protein_reference_to_node(protein_id)
            # print(protein_node)
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
                    # print(residue_node)
                    if residue_node is not None:
                        nodes.append(residue_node)

                        edges.append(
                            (residue_node[0], protein_id)
                        )
                        flag_node = self.data_.flag_to_node(residue, residue_node[0])
                        # print(flag_node)
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
                            flag_node = self.data_.flag_to_node(
                                flag, region_id)
                            if flag_node is not None:
                                nodes.append(flag_node)
                                edges.append(
                                    (flag_node[0], region_id)
                                )
                        for residue in region_data["residues"]:
                            residue_node = self.data_.residue_to_node(
                                residue, region_id)
                            if residue_node is not None:
                                nodes.append(residue_node)
                                edges.append(
                                    (residue_node[0], region_id)
                                )
                                edges.append(
                                    (residue_node[0], protein_id)
                                )

                                flag_node = self.data_.flag_to_node(
                                    residue, residue_node[0])
                                if flag_node is not None:
                                    nodes.append(flag_node)
                                    edges.append(
                                        (flag_node[0], residue_node[0])
                                    )
        graph.add_nodes_from(nodes)
        graph.add_edges_from(edges)

    def generate_families(self, families, graph):
        """Generate the nodes in the action graph for families."""
        print("Generating nodes for families...")
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
        """Generate the nodes in the action graph for small molecules."""
        print("Generating nodes for small molecules...")
        nodes = []
        for molecule_id in small_molecules:
            molecule_node = self.data_.small_molecule_to_node(molecule_id)
            if molecule_node is not None:
                nodes.append(molecule_node)
        graph.add_nodes_from(nodes)

    def generate_complexes(self, complexes, graph):
        """Generate the nodes in the action graph for complexes."""
        print("Generating nodes for complexes...")
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

    def generate_modification(self, modifications, graph):
        """Generate the nodes in the action graph and the nuggets for modifications."""
        print("Generating nodes for modifications...")

        nuggets = []

        for reaction_id, data in modifications.items():
            nugget = TypedDiGraph()

            # add source of modification the nugget
            # in addition get features of the sources
            # (reaction preconditions)
            for entity, features in data["sources"].items():
                if entity not in nugget.nodes():
                    nugget.add_node(
                        entity,
                        graph.node[entity].type_,
                        graph.node[entity].attrs_)
                for feature in features:
                    if self.data_.is_flag(feature):
                        flag_node = self.data_.flag_to_node(feature, entity)
                        modification = [k for k in flag_node[2].keys()][0]
                        nugget.add_node(
                            flag_node[0],
                            flag_node[1],
                            flag_node[2]
                        )
                        nugget.node[flag_node[0]].attrs_[modification] = 1
                        nugget.add_edge(flag_node[0], entity)
                    elif self.data_.is_residue(feature):
                        residue_node = self.data_.residue_to_node(feature, entity)
                        nugget.add_node(
                            residue_node[0],
                            residue_node[1],
                            residue_node[2]
                        )
                        nugget.add_edge(residue_node[0], entity)
                        flag_node = self.data_.flag_to_node(feature, residue_node[0])
                        modification = [k for k in flag_node[2].keys()][0]
                        nugget.add_node(
                            flag_node[0],
                            flag_node[1],
                            flag_node[2]
                        )
                        nugget.node[flag_node[0]].attrs_[modification] = 1
                        nugget.add_edge(flag_node[0], residue_node[0])
                    elif self.data_.is_fragment_feature(feature):
                        continue
                    else:
                        raise ValueError(
                            "Invalid modification %s!" % reaction_id)

            for entity, flag, value in data["targets"]:
                # add modification action to the
                # action graph
                target = str(entity) + "_state_" + str(flag)
                mod_node = self.data_.modification_to_node(reaction_id, target)
                mod_s = generate_modification_source(reaction_id, target)
                mod_t = generate_modification_target(reaction_id, target)

                graph.add_node(mod_node[0], mod_node[1])
                graph.add_node(mod_s[0], mod_s[1])
                graph.add_node(mod_t[0], mod_t[1])

                graph.add_edge(mod_s[0], mod_node[0])
                graph.add_edge(mod_t[0], mod_node[0])

                # add modification action to the nugget
                nugget.add_nodes_from([mod_node, mod_s, mod_t])

                nugget.add_edges_from([
                    (mod_s[0], mod_node[0]),
                    (mod_t[0], mod_node[0])
                ])

                # connect sources of modification with action
                # in the action graph
                for source in data["sources"].keys():
                    if source in graph.nodes():
                        graph.add_edge(source, mod_s[0])
                        nugget.add_edge(source, mod_s[0])
                    else:
                        print(source)

                if self.data_.is_flag(flag):
                    flag_node = self.data_.flag_to_node(flag, entity)
                    if entity not in nugget.nodes():
                        protein_node = self.data_.protein_reference_to_node(entity)
                        nugget.add_node(protein_node[0], protein_node[1], protein_node[2])
                    if flag_node[0] not in graph.nodes():
                        graph.add_node(flag_node[0], flag_node[1], flag_node[2])
                    graph.add_edge(mod_t[0], flag_node[0])
                    modification = [k for k in flag_node[2].keys()][0]
                    nugget.add_node(
                        flag_node[0],
                        flag_node[1],
                        flag_node[2])
                    nugget.node[flag_node[0]].attrs_[modification] = value
                    nugget.node[mod_node[0]].attrs_[modification] = value
                    nugget.add_edge(flag_node[0], protein_node[0])
                    nugget.add_edge(
                        mod_t[0], flag_node[0])
                elif self.data_.is_residue(flag):
                    residue_node = self.data_.residue_to_node(flag, entity)
                    nugget.add_node(residue_node[0], residue_node[1], residue_node[2])
                    if entity not in nugget.nodes():
                        protein_node = self.data_.protein_reference_to_node(entity)
                        nugget.add_node(protein_node[0], protein_node[1], protein_node[2])
                    nugget.add_edge(residue_node[0], protein_node[0])
                    flag_node = self.data_.flag_to_node(flag, residue_node[0])
                    if flag_node[0] not in graph.nodes():
                        graph.add_node(flag_node[0], flag_node[1], flag_node[2])
                    graph.add_edge(mod_t[0], flag_node[0])
                    modification = [k for k in flag_node[2].keys()][0]
                    nugget.add_node(
                        flag_node[0],
                        flag_node[1],
                        flag_node[2])
                    nugget.node[flag_node[0]].attrs_[modification] = value
                    nugget.node[mod_node[0]].attrs_[modification] = value
                    nugget.add_edge(flag_node[0], residue_node[0])
                    nugget.add_edge(
                        mod_t[0], flag_node[0])

                elif self.data_.is_fragment_feature(flag):
                    continue
                else:
                    raise ValueError(
                        "Invalid modification %s!" % reaction_id)
            nuggets.append(nugget)
        return nuggets

    def collect_agents(self, ignore_families=False):
        """Collect all the agents for the action graph construction."""
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
        """Collect all the interactions for the action graph and nuggets."""
        modifications = self.collect_modifications(ignore_families)
        # bindings = self.collect_bindings(ignore_families)
        actions = {
            "MOD": modifications,
            # "BND": bindings
        }
        return actions

    def import_model(self, filename, ignore_families=False):
        """Collect the data from BioPAX and generate action graph and nuggets."""
        self.data_.load(filename)
        graph = TypedDiGraph(self.metamodel_)

        agents = self.collect_agents(ignore_families)
        actions = self.collect_actions(ignore_families)
        self.generate_proteins(agents["proteins"], graph)
        self.generate_families(agents["protein_families"], graph)
        self.generate_small_molecules(agents["small_molecules"], graph)
        self.generate_families(agents["small_molecule_families"], graph)
        self.generate_complexes(agents["complexes"], graph)
        self.generate_families(agents["complex_families"], graph)
        modification_nuggets = self.generate_modification(
            actions["MOD"], graph)

        return graph, modification_nuggets
