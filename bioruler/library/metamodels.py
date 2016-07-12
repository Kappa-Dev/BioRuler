from regraph.library.data_structures import TypedDiGraph


metametamodel_AG = TypedDiGraph()

metametamodel_AG.add_nodes_from(
    [
        ("agent", "node"),
        ("action", "node"),
        ("attribute_node", "node")
    ]
)
metametamodel_AG.add_edges_from(
    [
        ("agent", "agent"),
        ("agent", "action"),
        ("action", "action"),
        ("action", "agent"),
        ("action", "attribute_node"),
        ("attribute_node", "attribute_node"),
        ("attribute_node", "agent"),
    ]
)

metamodel_AG = TypedDiGraph(metametamodel_AG)

metamodel_AG.add_nodes_from(
    [
        ("protein", "agent"),
        ("region", "agent"),
        ("residue", "agent"),
        ("family", "agent"),
        ("complex", "agent"),
        ("small_molecule", "agent"),
        ("flag", "attribute_node"),

        ("FAM", "action"),
        ("FAM_s", "action"),
        ("FAM_t", "action"),

        ("COM", "action"),
        ("COM_s", "action"),
        ("COM_t", "action"),

        ("BND", "action"),
        ("BND_s", "action"),

        ("BRK", "action"),
        ("BRK_t", "action"),

        ("MOD", "action"),
        ("MOD_s", "action"),
        ("MOD_t", "action")
    ])

metamodel_AG.add_edges_from(
    [
        ("region", "protein"),
        ("residue", "protein"),
        ("residue", "region"),
        ("flag", "protein"),
        ("flag", "region"),
        ("flag", "residue"),
        ("flag", "family"),
        ("flag", "complex"),
        ("flag", "small_molecule"),

        ("BND_s", "BND"),
        ("protein", "BND_s"),
        ("region", "BND_s"),
        ("small_molecule", "BND_s"),
        ("family", "BND_s"),
        ("complex", "BND_s"),

        ("BRK", "BND"),
        ("BRK_t", "BRK"),
        ("BRK_t", "protein"),
        ("BRK_t", "region"),
        ("BRK_t", "family"),
        ("BRK_t", "complex"),
        ("BRK_t", "small_molecule"),

        ("MOD_s", "MOD"),
        ("MOD_t", "MOD"),
        ("protein", "MOD_s"),
        ("region", "MOD_s"),
        ("family", "MOD_s"),
        ("complex", "MOD_s"),
        ("small_molecule", "MOD_s"),
        ("MOD_t", "flag"),

        ("FAM_s", "FAM"),
        ("FAM_t", "FAM"),
        ("protein", "FAM_s"),
        ("region", "FAM_s"),
        ("small_molecule", "FAM_s"),
        ("family", "FAM_s"),
        ("complex", "FAM_s"),
        ("FAM_t", "family"),

        ("COM_s", "COM"),
        ("COM_t", "COM"),
        ("protein", "COM_s"),
        ("small_molecule", "COM_s"),
        ("family", "COM_s"),
        ("complex", "COM_s"),
        ("COM_t", "complex")
    ]
)

simple_metamodel_AG = TypedDiGraph(metametamodel_AG)

simple_metamodel_AG.add_nodes_from(
    [
        ("protein", "agent"),
        ("region", "agent"),
        ("residue", "agent"),
        ("family", "agent"),
        ("complex", "agent"),
        ("small_molecule", "agent"),
        ("flag", "attribute_node"),
        ("FAM", "action"),
        ("COM", "action"),
        ("BND", "action"),
        ("BRK", "action"),
        ("MOD", "action"),
    ]
)

simple_metamodel_AG.add_edges_from(
    [
        ("region", "protein"),
        ("residue", "protein"),
        ("residue", "region"),
        ("flag", "protein"),
        ("flag", "region"),
        ("flag", "residue"),
        ("flag", "family"),
        ("flag", "complex"),
        ("flag", "small_molecule"),

        ("protein", "BND"),
        ("region", "BND"),
        ("small_molecule", "BND"),
        ("family", "BND"),
        ("complex", "BND"),

        ("BRK", "BND"),
        ("BRK", "protein"),
        ("BRK", "region"),
        ("BRK", "family"),
        ("BRK", "complex"),
        ("BRK", "small_molecule"),

        ("protein", "MOD"),
        ("region", "MOD"),
        ("family", "MOD"),
        ("complex", "MOD"),
        ("small_molecule", "MOD"),
        ("MOD", "flag"),

        ("protein", "FAM"),
        ("region", "FAM"),
        ("small_molecule", "FAM"),
        ("family", "FAM"),
        ("complex", "FAM"),
        ("FAM", "family"),

        ("protein", "COM"),
        ("small_molecule", "COM"),
        ("family", "COM"),
        ("complex", "COM"),
        ("COM", "complex")
    ]
)
