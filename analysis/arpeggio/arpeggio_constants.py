# Refer to https://github.com/harryjubb/arpeggio for more details.

CONTACTS_FIELDS = [
    "atom_1",
    "atom_2",
    "clash",
    "covalent",
    "vdw_clash",
    "vdw",
    "proximal",
    "hydrogen_bond",
    "weak_hydrogen_bond",
    "halogen_bond",
    "ionic",
    "metal_complex",
    "aromatic",
    "hydrophobic",
    "carbonyl",
    "polar",
    "weak_polar",
    "interacting_entities",
]

CHOSEN_CLASHES = [
    "clash",
    "vdw_clash",
]

CHOSEN_CONTACT_TYPES = [
    # "covalent",
    "vdw",
    # "proximal",
    "hydrogen_bond",
    "weak_hydrogen_bond",
    "halogen_bond",
    "ionic",
    "metal_complex",
    "aromatic",
    "hydrophobic",
    "carbonyl",
    "polar",
    "weak_polar",
]

CHOSEN_ARI_TYPES = [
    "CARBONPI",
    "CATIONPI",
    "DONORPI",
    "HALOGENPI",
    "METSULPHURPI"
]

assert set(CHOSEN_CONTACT_TYPES+CHOSEN_CLASHES).issubset(set(CONTACTS_FIELDS)), (
    "CHOSEN_CONTACT_TYPES and CHOSEN_CLASHES must be a subset of CONTACTS_FIELDS."
)

RI_FIELDS = [
    "ring1_id",
    "ring1_residue",
    "ring1_centroid",
    "ring2_id",
    "ring2_residue",
    "ring2_centroid",
    "interaction_type",
    "interacting_entites",
]

ARI_FIELDS = [
    "atom",
    "ring_id",
    "residue",
    "ring_centroid",
    "interaction_type",
]

PI_PI_STACKING_TYPES = [
    "FF",
    "OF",
    "EE",
    "FT",
    "OT",
    "ET",
    "FE",
    "OE",
    "EF",
]

RING_FIELDS = [
    "ring_id",
    "ring_residue",
    "ring_centroid",
]

INTERACTION_DF_COLS = {
    "atom": [
        "chain_1",
        "res_1",
        "chain_2",
        "res_2",
        "atom_1",
        "atom_2",
        "interaction_type",
    ],
    "residue": [
        "chain_1",
        "res_1",
        "chain_2",
        "res_2",
        "interaction_type",
    ],
}

MARKER_COMMAND = "marker #$marker_model_id position $x,$y,$z color $marker_color radius $marker_radius"
PBOND_COMMAND = "pbond #$spec1#$spec2 color $bond_color radius $bond_radius name $bond_name dashes $bond_dashes"
TRANSPARENCY_COMMAND = "transparency #$model_id $transparency target $target_spec"

MARKER_ATTRIBUTES = {
    "model_id": 2,
    "color": "white",
    "radius": 0.1,
}

PBOND_ATTRIBUTES = {
    "default": {
        "color": "yellow",
        "radius": 0.1,
        "dashes": 6,
        "transparency": 0,
    },
    "vdw_clash": {
        "color": "magenta",
        "radius": 0.2,
        "dashes": 0,
        "transparency": 0,
    },
    "clash": {
        "color": "magenta",
        "radius": 0.2,
        "dashes": 0,
        "transparency": 0,
    },
    "covalent": {
        "color": "gray",
        "radius": 0.2,
        "dashes": 0,
        "transparency": 0,
    },
    "vdw": {
        "color": "lightgray",
        "radius": 0.1,
        "dashes": 6,
        "transparency": 0,
    },
    "proximal": {
        "color": "white",
        "radius": 0.1,
        "dashes": 6,
        "transparency": 90,
    },
    "hydrogen_bond": {
        "color": "blue",
        "radius": 0.1,
        "dashes": 10,
        "transparency": 0,
    },
    "weak_hydrogen_bond": {
        "color": "blue",
        "radius": 0.1,
        "dashes": 10,
        "transparency": 80,
    },
    "halogen_bond": {
        "color": "cyan",
        "radius": 0.1,
        "dashes": 10,
        "transparency": 0,
    },
    "ionic": {
        "color": "red",
        "radius": 0.1,
        "dashes": 6,
        "transparency": 0,
    },
    "metal_complex": {
        "color": "hotpink",
        "radius": 0.1,
        "dashes": 30,
        "transparency": 0,
    },
    "aromatic": {
        "color": "green",
        "radius": 0.15,
        "dashes": 20,
        "transparency": 0,
    },
    "hydrophobic": {
        "color": "orange",
        "radius": 0.1,
        "dashes": 6,
        "transparency": 0,
    },
    "carbonyl": {
        "color": "purple",
        "radius": 0.1,
        "dashes": 6,
        "transparency": 0,
    },
    "polar": {
        "color": "cyan",
        "radius": 0.1,
        "dashes": 6,
        "transparency": 0,
    },
    "weak_polar": {
        "color": "cyan",
        "radius": 0.1,
        "dashes": 6,
        "transparency": 80,
    },
    "ri": {
        "color": "lime",
        "radius": 0.3,
        "dashes": 30,
        "transparency": 0,
    },
    "ari": {
        "color": "cornflowerblue",
        "radius": 0.1,
        "dashes": 30,
        "transparency": 0,
    },
}

DOCKER_BASE_COMMAND = (
    """sudo docker run --rm -v $path_to_mount:$container_path -u `id -u`:`id -g` -it harryjubb/arpeggio python arpeggio.py $input_pdb_path -v -wh
    """
)

MODEL_ID = 1