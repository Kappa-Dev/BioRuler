import warnings

from jpype import (java, startJVM, getDefaultJVMPath, JPackage)


def generate_family_action_node(family_id):
    """."""
    node = {
        "id": "%s_family" % family_id,
        "type": "FAM",
    }
    return node


def generate_family_source(family_id):
    """."""
    node = {
        "id": "%s_family_s" % family_id,
        "type": "FAM_s",
    }
    return node


def generate_family_target(family_id):
    """."""
    node = {
        "id": "%s_family_t" % family_id,
        "type": "FAM_t",
    }
    return node


def edge_from_ids(source, target, attributes=None):
    """."""
    edge = {
        "from": source,
        "to": target
    }
    if attributes is not None:
        edge["attrs"] = attributes
    return edge


class BioPAXModel():

    def __init__(self):
        """Initialize Java and load BioPAX classes of interest."""
        startJVM(getDefaultJVMPath(), "-ea", "-Xmx1g", "-Djava.class.path=paxtools-5.0.0-20160601.jar")
        self.java_io_ = JPackage("java.io")
        self.paxtools_ = JPackage("org.biopax.paxtools")
        self.io_ = self.paxtools_.io.SimpleIOHandler(self.paxtools_.model.BioPAXLevel.L3)

        # That's how we start extracting agents
        # Protein
        # Small molecule
        # Rna
        # Complex
        # Dna
        self.protein_class_ = java.lang.Class.forName(
            "org.biopax.paxtools.model.level3.Protein", True,
            java.lang.ClassLoader.getSystemClassLoader())
        self.protein_reference_class_ = java.lang.Class.forName(
            "org.biopax.paxtools.model.level3.ProteinReference", True,
            java.lang.ClassLoader.getSystemClassLoader())
        self.fragment_feature_class_ = java.lang.Class.forName(
            "org.biopax.paxtools.model.level3.FragmentFeature", True,
            java.lang.ClassLoader.getSystemClassLoader())
        self.modification_feature_class_ = java.lang.Class.forName(
            "org.biopax.paxtools.model.level3.ModificationFeature", True,
            java.lang.ClassLoader.getSystemClassLoader())
        self.small_molecule_class_ = java.lang.Class.forName(
            "org.biopax.paxtools.model.level3.SmallMolecule", True,
            java.lang.ClassLoader.getSystemClassLoader())
        self.small_molecule_reference_class_ = java.lang.Class.forName(
            "org.biopax.paxtools.model.level3.SmallMoleculeReference", True,
            java.lang.ClassLoader.getSystemClassLoader())
        self.rna_class_ = java.lang.Class.forName(
            "org.biopax.paxtools.model.level3.Rna", True,
            java.lang.ClassLoader.getSystemClassLoader())
        self.rna_reference_class_ = java.lang.Class.forName(
            "org.biopax.paxtools.model.level3.RnaReference", True,
            java.lang.ClassLoader.getSystemClassLoader())
        self.complex_class_ = java.lang.Class.forName(
            "org.biopax.paxtools.model.level3.Complex", True,
            java.lang.ClassLoader.getSystemClassLoader())
        self.dna_class_ = java.lang.Class.forName(
            "org.biopax.paxtools.model.level3.Dna", True,
            java.lang.ClassLoader.getSystemClassLoader())
        self.dna_reference_class_ = java.lang.Class.forName(
            "org.biopax.paxtools.model.level3.DnaReference", True,
            java.lang.ClassLoader.getSystemClassLoader())
        self.model_ = None

    def load(self, filename):
        """Import a BioPAX model from the file."""
        file_is = self.java_io_.FileInputStream(
            filename)
        self.model_ = self.io_.convertFromOWL(file_is)
        file_is.close()

    def is_protein_family(self, reference_id):
        """."""
        reference = self.model_.getByID(reference_id)
        if len(reference.getMemberEntityReference()) > 0:
            return True
        else:
            physical_entities = reference.getEntityReferenceOf()
            for entity in physical_entities:
                if len(entity.getMemberPhysicalEntity()) > 0:
                    return True
            return False

    def is_complex_family(self, complex_id):
        """."""
        complex = self.model_.getByID(complex_id)
        return len(complex.getMemberPhysicalEntity()) > 0

    def is_fragment(self, protein_id):
        """."""
        protein = self.model_.getByID(protein_id)
        features = protein.getFeature()
        for f in features:
            if f.getModelInterface() == self.fragment_feature_class_:
                return True
        return False

    def protein_reference_to_node(self, protein_reference_id):
        """."""

        protein_reference = self.model_.getByID(protein_reference_id)
        xref = set(protein_reference.getXref())
        if len(xref) > 1:
            warnings.warn(
                "Protein reference %s (%s) has ambiguous Unified Reference!" %
                (str(protein_reference.getName()), str(protein_reference)))
        elif len(xref) < 1:
            warnings.warn(
                "Protein reference %s (%s) does not have Unified Reference!" %
                (str(protein_reference.getName()), str(protein_reference)))

        protein_attrs = {}
        if len(xref) == 1:
            protein_id = list(xref)[0]
            protein_attrs[protein_id.getDb()] = protein_id.getId()
        if len(xref) > 1:
            for el in xref:
                protein_attrs[el.getDb()] = el.getId()

        protein_attrs["Name"] = list(protein_reference.getName())
        locations = set()
        for entity in protein_reference.getEntityReferenceOf():
            location = entity.getCellularLocation()
            if location is not None:
                locations.update(list(location.getTerm()))
        if len(locations) > 0:
            protein_attrs["loc"] = list(locations)
        node = {
            "id": protein_reference.getUri(),
            "type": "protein",
            "attrs": protein_attrs
        }
        return node

    def region_to_node(self, region_feature_id):
        """."""
        region_feature = self.model_.getByID(region_feature_id)
        location = region_feature.getFeatureLocation()
        region_id = region_feature.getUri()
        region_attrs = {}
        start = location.getSequenceIntervalBegin()
        if start is not None:
            region_attrs["start"] = start.getSequencePosition()
        end = location.getSequenceIntervalEnd()
        if end is not None:
            region_attrs["end"] = end.getSequencePosition()

        node = {
            "id": region_id,
            "type": "region",
            "attrs": region_attrs
        }
        return node

    def residue_to_node(self, protein_id, residue_id):
        """."""
        residue = self.model_.getByID(residue_id)
        residue_attrs = {}
        residue_attrs["loc"] = residue.getFeatureLocation().getSequencePosition()
        node = {
            "id": "%s@%s" % (protein_id, residue_attrs["loc"]),
            "type": "residue",
            "attrs": residue_attrs,
        }
        return node

    def flag_to_node(self, entitiy_id, flag_id):
        """."""
        flag = self.model_.getByID(flag_id)
        flag_attrs = {}
        states = list(flag.getModificationType().getTerm())
        if len(states) == 1:
            flag_attrs[states[0]] = [0, 1]
            node = {
                "id": "%s_%s_flag" % (entitiy_id, flag_id),
                "type": "flag",
                "attrs": flag_attrs,
            }
            return node
        else:
            warnings.warn("Ambiguous state (%s)! Cannot convert to node" % states)

    def family_reference_to_node(self, family_reference_id):
        """."""
        family_reference = self.model_.getByID(family_reference_id)
        family_attrs = {}
        family_attrs["Name"] = list(family_reference.getName())
        node = {
            "id": family_reference.getUri(),
            "type": "family",
            "attrs": family_attrs
        }
        return node

    def small_molecule_to_node(self, small_molecule_id):
        """."""
        small_molecule_reference = self.model_.getByID(small_molecule_id)
        xref = set(small_molecule_reference.getXref())
        if len(xref) > 1:
            warnings.warn(
                "Small molecule reference %s (%s) has ambiguous Unified Reference!" %
                (str(small_molecule_reference.getName()), str(small_molecule_reference)))
        elif len(xref) == 0:
            warnings.warn(
                "Small molecule reference %s (%s) does not have Unified Reference!" %
                (str(small_molecule_reference.getName()), str(small_molecule_reference)))

        molecule_attrs = {}
        if len(xref) == 1:
            molecule_id = list(xref)[0]
            molecule_attrs[molecule_id.getDb()] = molecule_id.getId()
        elif len(xref) > 1:
            for el in xref:
                molecule_attrs[el.getDb()] = el.getId()
        molecule_attrs["Name"] = list(small_molecule_reference.getName())
        physical_entities = small_molecule_reference.getEntityReferenceOf()
        locations = set()
        for entity in physical_entities:
            location = entity.getCellularLocation()
            if location is not None:
                locations.update(list(location.getTerm()))
        if len(locations) > 0:
            molecule_attrs["loc"] = list(locations)
        node = {
            "id": small_molecule_reference.getUri(),
            "type": "small_molecule",
            "attrs": molecule_attrs
        }
        return node

    def complex_to_node(self, complex_id):
        complex = self.model_.getByID(complex_id)
        complex_attrs = {}
        complex_attrs["Name"] = list(complex.getName())
        if complex.getCellularLocation() is not None:
            complex_attrs["loc"] = list(complex.getCellularLocation().getTerm())
        node = {
            "id": complex.getUri(),
            "type": "complex",
            "attrs": complex_attrs
        }
        return node

    def get_modifications(self, physical_entities):
        """."""
        residues = set()
        flags = set()
        for e in physical_entities:
            entity = self.model_.getByID(e)
            features = entity.feature
            for f in features:
                if f.getModelInterface() == self.modification_feature_class_:
                    if f.getFeatureLocation() is not None:
                        residues.add(f.getUri())
                    else:
                        flags.add(f.getUri())
        return {"residues": residues, "flags": flags}

    def residue_in_region(self, residue_id, region_id):
        """."""
        residue = self.model_.getByID(residue_id)
        region = self.model_.getByID(region_id)

        region_location = region.getFeatureLocation()
        region_start = region_location.getSequenceIntervalBegin()
        region_end = region_location.getSequenceIntervalEnd()

        residue_location = residue.getFeatureLocation()
        if region_start is not None and region_end is not None and residue_location is not None:
            x = residue_location.getSequencePosition()
            start = region_start.getSequencePosition()
            end = region_end.getSequencePosition()
            return (x >= start) and (x <= end)
        else:
            return False
