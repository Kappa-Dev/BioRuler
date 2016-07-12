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
        ("protein", "complex"),
        ("residue", "complex"),
        ("small_molecule", "complex"),
        ("family", "complex"),
        ("complex", "complex"),

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
    ]
)