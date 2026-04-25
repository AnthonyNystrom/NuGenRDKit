"""
Phase-0 inventory — the single source of truth for regression tests.

Every HTML page, every /api/v1/** endpoint, every uploadable sample file, and
every "Load Example" button in the app is listed here. Tests consume this
module; if a route is added/removed without a matching inventory entry, the
coverage test fails.

This file is the oracle for the rest of the upgrade: Phase 1 (backend
hardening) and Phase 2 (code-quality refactor) must not change any payload
shape listed here, and Phase 5 (page migration) must not change any page URL.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


def _build_simple_rxn_block() -> str:
    """
    Build a minimal but valid MDL RXN block for upload_rxn tests.
    RDKit round-trip: SMARTS → Reaction → RxnBlock ensures formatting stays
    correct across RDKit versions.
    """
    from rdkit.Chem import AllChem

    rxn = AllChem.ReactionFromSmarts(
        "[C:1](=[O:2])[OH:3].[N:4]>>[C:1](=[O:2])[N:4].[O:3]",
        useSmiles=False,
    )
    return AllChem.ReactionToRxnBlock(rxn)


# --------------------------------------------------------------------------- #
# Canonical sample molecules
# --------------------------------------------------------------------------- #

ASPIRIN_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"
CAFFEINE_SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
ETHANOL_SMILES = "CCO"
BENZENE_SMILES = "c1ccccc1"
REACTION_SMARTS = "[C:1](=[O:2])[OH:3].[N:4]>>[C:1](=[O:2])[N:4].[O:3]"
ATTACK_SMARTS = "[OH]"  # substructure pattern for phenol

ETHANOL_MOLBLOCK = """
  Mrv2311 01012023000000000000

  3  2  0  0  0  0            999 V2000
   -0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
"""

# A minimal valid InChI for ethanol.
ETHANOL_INCHI = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"

# A one-molecule SDF (ethanol) — used for /structure/sdf_to_smiles etc.
ETHANOL_SDF = ETHANOL_MOLBLOCK.lstrip("\n") + "$$$$\n"

# A minimal PDB block for methane (one atom is fine for round-trip tests).
METHANE_PDB = (
    "HETATM    1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n"
)

# .smi content for file-upload tests (newline-separated SMILES + optional name).
SMI_FILE_BODY = f"{ASPIRIN_SMILES} aspirin\n{CAFFEINE_SMILES} caffeine\n{ETHANOL_SMILES} ethanol\n"

# A minimal valid .sdf body (single molecule, ethanol) for file-import tests.
SDF_FILE_BODY = ETHANOL_MOLBLOCK.lstrip("\n") + "$$$$\n"

# A minimal valid .rxn body generated at module import time. Kept in sync
# with whatever RDKit version is installed, so upload_rxn tests don't drift.
RXN_FILE_BODY = _build_simple_rxn_block()


# --------------------------------------------------------------------------- #
# Pages (HTML routes)
# --------------------------------------------------------------------------- #

PAGES: list[str] = [
    "/",
    "/structure",
    "/descriptors",
    "/fingerprints",
    "/similarity",
    "/coordinates",
    "/properties",
    "/reactions",
    "/visualization",
]

# Non-HTML GET routes that should still respond 200 with `success: True`.
JSON_GET_ROUTES: list[str] = [
    "/api",
    "/health",
    "/api/v1/descriptors/list",
    # Phase 7: diagnostic probes (mounted under /api/status, not /api/v1).
    "/api/status/rdkit",
    "/api/status/disk",
    "/api/status/",
]


# Phase-7 page-style routes (return HTML). Tested by test_routes_ok.
HTML_DOC_ROUTES: list[str] = [
    "/api/v1/docs",  # Swagger UI
]


# Phase-7 GET routes that return JSON but DON'T follow the standard
# `success: True` envelope (e.g. the OpenAPI spec is its own format).
OTHER_JSON_GET_ROUTES: list[str] = [
    "/api/v1/openapi.json",
]


# --------------------------------------------------------------------------- #
# API endpoints — POST payloads known to succeed on the current app
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Endpoint:
    """One /api/v1/** smoke-test case."""

    url: str
    payload: dict[str, Any] = field(default_factory=dict)
    # If True this URL is an alias for another handler and only needs to
    # resolve (we don't assert on response shape to avoid double-pinning).
    alias: bool = False
    # Keys that must be present in the JSON response on success. Defaults
    # to {"success"}. Phase 1 will standardize this further.
    expect_keys: tuple[str, ...] = ("success",)
    # Set to True when the endpoint is too slow / too complex for the
    # Phase-0 smoke suite and should only run under pytest --slow.
    slow: bool = False
    # Skip the smoke assertion but still register the URL for coverage.
    # Used for endpoints we know are broken today; tests document them.
    known_broken: str | None = None


# --- /api/v1/structure ---------------------------------------------------- #

STRUCTURE_ENDPOINTS: list[Endpoint] = [
    Endpoint(
        "/api/v1/structure/smiles_to_mol",
        {"smiles": ETHANOL_SMILES},
        expect_keys=("success", "mol_block"),
    ),
    Endpoint("/api/v1/structure/mol", {"smiles": ETHANOL_SMILES}, alias=True),
    Endpoint(
        "/api/v1/structure/mol_to_smiles",
        {"mol_block": ETHANOL_MOLBLOCK},
        expect_keys=("success", "smiles", "canonical_smiles"),
    ),
    Endpoint(
        "/api/v1/structure/inchi_to_smiles",
        {"inchi": ETHANOL_INCHI},
        expect_keys=("success", "smiles"),
    ),
    Endpoint(
        "/api/v1/structure/smiles_to_inchi",
        {"smiles": ETHANOL_SMILES},
        expect_keys=("success", "inchi", "inchi_key"),
    ),
    Endpoint("/api/v1/structure/inchi", {"smiles": ETHANOL_SMILES}, alias=True),
    Endpoint(
        "/api/v1/structure/canonicalize",
        {"smiles": "OCC"},
        expect_keys=("success", "canonical_smiles"),
    ),
    Endpoint(
        "/api/v1/structure/validate",
        {"smiles": ETHANOL_SMILES},
        expect_keys=("success", "is_valid"),
    ),
    Endpoint("/api/v1/structure/standardize", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/structure/add_hydrogens", {"smiles": ETHANOL_SMILES}),
    Endpoint("/api/v1/structure/remove_hydrogens", {"smiles": ETHANOL_SMILES}),
    Endpoint("/api/v1/structure/stereochemistry", {"smiles": "C[C@H](N)C(=O)O"}),
    Endpoint(
        "/api/v1/structure/enumerate_stereoisomers", {"smiles": "C[C@H](N)C(=O)O", "max_isomers": 4}
    ),
    Endpoint("/api/v1/structure/kekulize", {"smiles": BENZENE_SMILES}),
    Endpoint("/api/v1/structure/aromaticity", {"smiles": BENZENE_SMILES}),
    Endpoint("/api/v1/structure/fragment", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/structure/murcko_scaffold", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/structure/ring_info", {"smiles": BENZENE_SMILES}),
    # /sdf_to_smiles calls MolFromMolBlock under the hood, so it wants a
    # plain mol block — not a $$$$-terminated SDF. Phase 2 may rename or
    # fix the handler; for now we lock the current contract.
    Endpoint("/api/v1/structure/sdf_to_smiles", {"sdf_block": ETHANOL_MOLBLOCK}),
    Endpoint("/api/v1/structure/smiles_to_sdf", {"smiles": ETHANOL_SMILES}),
    Endpoint("/api/v1/structure/pdb_to_smiles", {"pdb_block": METHANE_PDB}),
    Endpoint("/api/v1/structure/smiles_to_pdb", {"smiles": ETHANOL_SMILES}),
    Endpoint("/api/v1/structure/neutralize", {"smiles": "CC(=O)[O-]"}),
    Endpoint("/api/v1/structure/remove_fragments", {"smiles": f"{ETHANOL_SMILES}.O"}),
    Endpoint("/api/v1/structure/cleanup", {"smiles": ETHANOL_SMILES}),
]

# --- /api/v1/descriptors -------------------------------------------------- #

DESCRIPTORS_ENDPOINTS: list[Endpoint] = [
    Endpoint("/api/v1/descriptors/basic", {"smiles": ETHANOL_SMILES}, expect_keys=("success",)),
    Endpoint("/api/v1/descriptors/lipinski", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/descriptors/logp", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/descriptors/topological", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/descriptors/all", {"smiles": ETHANOL_SMILES}, slow=True),
    Endpoint("/api/v1/descriptors/vsa", {"smiles": ASPIRIN_SMILES}),
]

# --- /api/v1/fingerprints ------------------------------------------------- #

FINGERPRINTS_ENDPOINTS: list[Endpoint] = [
    Endpoint(
        "/api/v1/fingerprints/morgan", {"smiles": ASPIRIN_SMILES, "radius": 2, "n_bits": 2048}
    ),
    Endpoint("/api/v1/fingerprints/rdkit", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/fingerprints/maccs", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/fingerprints/avalon", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/fingerprints/atom_pairs", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/fingerprints/topological_torsions", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/fingerprints/pattern", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/fingerprints/layered", {"smiles": ASPIRIN_SMILES}),
    Endpoint(
        "/api/v1/fingerprints/compare",
        {"smiles1": ETHANOL_SMILES, "smiles2": ASPIRIN_SMILES, "fp_type": "morgan"},
    ),
]

# --- /api/v1/similarity --------------------------------------------------- #

_SIM_PAIR = {"query_smiles": ETHANOL_SMILES, "target_smiles": ASPIRIN_SMILES}

SIMILARITY_ENDPOINTS: list[Endpoint] = [
    Endpoint("/api/v1/similarity/tanimoto", _SIM_PAIR),
    Endpoint("/api/v1/similarity/dice", _SIM_PAIR),
    Endpoint("/api/v1/similarity/cosine", _SIM_PAIR),
    Endpoint("/api/v1/similarity/sokal", _SIM_PAIR),
    Endpoint("/api/v1/similarity/braunblanquet", _SIM_PAIR),
    Endpoint("/api/v1/similarity/kulczynski", _SIM_PAIR),
    Endpoint("/api/v1/similarity/mcconnaughey", _SIM_PAIR),
    Endpoint("/api/v1/similarity/rogotgoldberg", _SIM_PAIR),
    Endpoint("/api/v1/similarity/russel", _SIM_PAIR),
    Endpoint(
        "/api/v1/similarity/tversky",
        {**_SIM_PAIR, "alpha": 0.5, "beta": 0.5},
    ),
    Endpoint("/api/v1/similarity/asymmetric", _SIM_PAIR),
    Endpoint("/api/v1/similarity/allbit", _SIM_PAIR),
    Endpoint("/api/v1/similarity/onbit", _SIM_PAIR),
    Endpoint(
        "/api/v1/similarity/bulk_similarity",
        {
            "query_smiles": ETHANOL_SMILES,
            "target_smiles_list": [ASPIRIN_SMILES, CAFFEINE_SMILES, BENZENE_SMILES],
            "fp_type": "morgan",
            "similarity_metric": "tanimoto",
            "threshold": 0.0,
        },
    ),
    Endpoint(
        "/api/v1/similarity/substructure_search",
        {
            "pattern_smiles": ATTACK_SMARTS,
            "target_smiles_list": [ASPIRIN_SMILES, CAFFEINE_SMILES, ETHANOL_SMILES],
        },
    ),
    Endpoint(
        "/api/v1/similarity/substructure",
        {
            "pattern_smiles": ATTACK_SMARTS,
            "target_smiles_list": [ASPIRIN_SMILES, CAFFEINE_SMILES, ETHANOL_SMILES],
        },
        alias=True,
    ),
    Endpoint(
        "/api/v1/similarity/maximum_common_substructure",
        {"smiles1": ASPIRIN_SMILES, "smiles2": "CC(=O)OC1=CC=CC=C1"},
    ),
    Endpoint(
        "/api/v1/similarity/diverse_subset",
        {
            "smiles_list": [
                ASPIRIN_SMILES,
                CAFFEINE_SMILES,
                ETHANOL_SMILES,
                BENZENE_SMILES,
                "CCN",
                "CCC",
                "CCCC",
                "CCCCC",
            ],
            "num_compounds": 3,
            "fp_type": "morgan",
        },
    ),
]

# --- /api/v1/coordinates -------------------------------------------------- #

COORDINATES_ENDPOINTS: list[Endpoint] = [
    Endpoint("/api/v1/coordinates/generate_2d", {"smiles": ETHANOL_SMILES}),
    Endpoint("/api/v1/coordinates/2d", {"smiles": ETHANOL_SMILES}, alias=True),
    Endpoint(
        "/api/v1/coordinates/generate_3d",
        {
            "smiles": ETHANOL_SMILES,
            "optimize": True,
            "num_conformers": 1,
            "force_field": "UFF",
            "max_iterations": 50,
        },
    ),
    Endpoint(
        "/api/v1/coordinates/3d",
        {"smiles": ETHANOL_SMILES, "num_conformers": 1, "force_field": "UFF"},
        alias=True,
    ),
    Endpoint(
        "/api/v1/coordinates/optimize_geometry",
        {"smiles": ETHANOL_SMILES, "force_field": "UFF", "max_iterations": 50},
    ),
    Endpoint(
        "/api/v1/coordinates/optimize",
        {"smiles": ETHANOL_SMILES, "force_field": "UFF"},
        alias=True,
    ),
    # align_molecules requires the probe to be a substructure of the
    # reference (or vice versa); aligning aspirin to itself is the simplest
    # payload that exercises the path.
    Endpoint(
        "/api/v1/coordinates/align_molecules",
        {"reference_smiles": ASPIRIN_SMILES, "probe_smiles": ASPIRIN_SMILES},
    ),
    Endpoint(
        "/api/v1/coordinates/o3a_align",
        {"reference_smiles": ASPIRIN_SMILES, "probe_smiles": ASPIRIN_SMILES},
    ),
    Endpoint(
        "/api/v1/coordinates/conformer_search",
        {"smiles": ETHANOL_SMILES, "num_conformers": 3, "optimize": False},
        slow=True,
    ),
    Endpoint(
        "/api/v1/coordinates/multiple",
        {"smiles": ETHANOL_SMILES, "num_conformers": 3},
        alias=True,
        slow=True,
    ),
    Endpoint(
        "/api/v1/coordinates/constrained_embed",
        {"smiles": ETHANOL_SMILES, "constraints": []},
        slow=True,
    ),
    Endpoint(
        "/api/v1/coordinates/geometry_transforms",
        {"smiles": ETHANOL_SMILES, "transform_type": "centroid"},
    ),
]

# --- /api/v1/properties --------------------------------------------------- #

PROPERTIES_ENDPOINTS: list[Endpoint] = [
    Endpoint("/api/v1/properties/physicochemical", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/properties/drug_likeness", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/properties/qed", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/properties/fragments", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/properties/brics", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/properties/recap", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/properties/scaffold", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/properties/aromaticity", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/properties/charges", {"smiles": ETHANOL_SMILES}),
    Endpoint("/api/v1/properties/stereochemistry", {"smiles": "C[C@H](N)C(=O)O"}),
    Endpoint("/api/v1/properties/admet", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/properties/all", {"smiles": ETHANOL_SMILES}, slow=True),
]

# --- /api/v1/reactions ---------------------------------------------------- #

REACTIONS_ENDPOINTS: list[Endpoint] = [
    # `/process` expects a flat list of reactant SMILES (one per reaction
    # template), not a list-of-lists. See process_smarts_reaction().
    Endpoint(
        "/api/v1/reactions/process",
        {
            "reaction_smarts": REACTION_SMARTS,
            "reactants": ["CC(=O)O", "CN"],
            "reaction_type": "smarts",
            "max_products": 5,
        },
    ),
    Endpoint("/api/v1/reactions/parse_smarts", {"reaction_smarts": REACTION_SMARTS}),
    Endpoint(
        "/api/v1/reactions/run_reaction",
        {
            "reaction_smarts": REACTION_SMARTS,
            "reactant_smiles": ["CC(=O)O", "CN"],
            "max_products": 5,
        },
    ),
    Endpoint(
        "/api/v1/reactions/validate_reaction",
        {
            "reactant_smiles": ["CC(=O)O", "CN"],
            "product_smiles": ["CC(=O)NC", "O"],
        },
    ),
    # reaction_center requires exactly one reactant + one product (see
    # routes/reactions.py:579). Use a benign self-mapping.
    Endpoint(
        "/api/v1/reactions/reaction_center",
        {
            "reactant_smiles": ["CC(=O)O"],
            "product_smiles": ["CC(=O)OC"],
        },
    ),
    Endpoint(
        "/api/v1/reactions/enumerate_library",
        {
            "reaction_smarts": REACTION_SMARTS,
            "reactant_lists": [["CC(=O)O"], ["CN", "CCN"]],
            "max_products": 10,
        },
    ),
    Endpoint(
        "/api/v1/reactions/reaction_fingerprint",
        {"reaction_smarts": REACTION_SMARTS, "fp_type": "structural"},
    ),
    Endpoint(
        "/api/v1/reactions/brics_react",
        {"molecules": [ASPIRIN_SMILES, CAFFEINE_SMILES], "max_products": 10},
    ),
    Endpoint("/api/v1/reactions/recap_react", {"smiles": ASPIRIN_SMILES}),
]

# --- /api/v1/visualization ------------------------------------------------ #

VISUALIZATION_ENDPOINTS: list[Endpoint] = [
    Endpoint("/api/v1/visualization/draw_svg", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/visualization/draw_png", {"smiles": ASPIRIN_SMILES}),
    Endpoint(
        "/api/v1/visualization/draw_grid",
        {
            "smiles_list": [ASPIRIN_SMILES, CAFFEINE_SMILES, ETHANOL_SMILES],
            "mols_per_row": 2,
            "format": "svg",
        },
    ),
    Endpoint(
        "/api/v1/visualization/draw_reaction",
        {"reaction_smarts": REACTION_SMARTS, "format": "svg"},
    ),
    Endpoint(
        "/api/v1/visualization/highlight_substructure",
        {"molecule_smiles": ASPIRIN_SMILES, "pattern_smarts": ATTACK_SMARTS, "format": "svg"},
    ),
    Endpoint("/api/v1/visualization/draw_3d", {"smiles": ETHANOL_SMILES}),
    Endpoint("/api/v1/visualization/3d", {"smiles": ETHANOL_SMILES}, alias=True),
    Endpoint("/api/v1/visualization/surface", {"smiles": ETHANOL_SMILES}, slow=True),
    Endpoint("/api/v1/visualization/pharmacophore", {"smiles": ASPIRIN_SMILES}),
    Endpoint("/api/v1/visualization/scaffold", {"smiles": ASPIRIN_SMILES}),
    # /similarity_map uses `smiles` + optional `ref_smiles` (not smiles1/2).
    Endpoint(
        "/api/v1/visualization/similarity_map",
        {
            "smiles": ASPIRIN_SMILES,
            "ref_smiles": CAFFEINE_SMILES,
            "fp_type": "morgan",
            "map_type": "similarity",
        },
    ),
    # fingerprint_bit requires a bit that is actually set in the fingerprint
    # for this molecule + parameters. 389 is one of the Morgan bits set for
    # aspirin at radius=2, n_bits=2048 (confirmed via the endpoint's own
    # "available_bits" suggestion list).
    Endpoint(
        "/api/v1/visualization/fingerprint_bit",
        {"smiles": ASPIRIN_SMILES, "bit_id": 389, "fp_type": "morgan", "radius": 2, "n_bits": 2048},
    ),
    Endpoint(
        "/api/v1/visualization/fingerprint_env",
        {"smiles": ASPIRIN_SMILES, "atom_index": 0, "radius": 2},
    ),
    # matrix_grid takes a 2-D SMILES array (rows × columns) and renders a
    # comparison grid.
    Endpoint(
        "/api/v1/visualization/matrix_grid",
        {"smiles_matrix": [[ETHANOL_SMILES, ASPIRIN_SMILES], [CAFFEINE_SMILES, BENZENE_SMILES]]},
    ),
]


ALL_ENDPOINTS: list[Endpoint] = [
    *STRUCTURE_ENDPOINTS,
    *DESCRIPTORS_ENDPOINTS,
    *FINGERPRINTS_ENDPOINTS,
    *SIMILARITY_ENDPOINTS,
    *COORDINATES_ENDPOINTS,
    *PROPERTIES_ENDPOINTS,
    *REACTIONS_ENDPOINTS,
    *VISUALIZATION_ENDPOINTS,
]


# --------------------------------------------------------------------------- #
# File-upload endpoints — multipart/form-data
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class UploadEndpoint:
    url: str
    filename: str = "molecules.smi"
    file_body: str = SMI_FILE_BODY
    form: dict[str, str] = field(default_factory=dict)
    file_field: str = "file"
    slow: bool = False
    known_broken: str | None = None


UPLOAD_ENDPOINTS: list[UploadEndpoint] = [
    UploadEndpoint(
        "/api/v1/structure/batch_convert",
        form={"output_format": "canonicalize"},
    ),
    UploadEndpoint(
        "/api/v1/descriptors/batch_file",
        form={"descriptor_type": "basic", "output_format": "json"},
    ),
    UploadEndpoint(
        "/api/v1/fingerprints/batch_file",
        form={"fingerprint_type": "morgan", "output_format": "json"},
    ),
    UploadEndpoint(
        "/api/v1/similarity/bulk_similarity_file",
        form={
            "query_smiles": ETHANOL_SMILES,
            "fp_type": "morgan",
            "similarity_metric": "tanimoto",
            "threshold": "0.0",
        },
        slow=True,
    ),
    UploadEndpoint(
        "/api/v1/similarity/substructure_search_file",
        form={"pattern_smiles": ATTACK_SMARTS},
    ),
    UploadEndpoint(
        "/api/v1/similarity/diverse_subset_file",
        form={"num_compounds": "2", "fp_type": "morgan"},
        slow=True,
    ),
    UploadEndpoint(
        "/api/v1/properties/batch_file",
        form={"property_type": "physicochemical", "output_format": "json"},
        slow=True,
    ),
    UploadEndpoint(
        "/api/v1/visualization/batch_file",
        form={"viz_type": "2d", "output_format": "grid"},
    ),
    UploadEndpoint(
        "/api/v1/coordinates/import_file",
        filename="molecule.sdf",
        file_body=SDF_FILE_BODY,
    ),
    UploadEndpoint(
        "/api/v1/reactions/process_file",
        form={
            "reaction_type": "enumeration",
            # Phenyl ring with attachment point — valid SMILES.
            "core_structure": "c1ccc([*:1])cc1",
            "max_products": "5",
        },
    ),
    # Use a single-reactant reaction (chloride hydrolysis) so the single-
    # file UploadEndpoint model covers it. Multi-file variants can be
    # exercised via a dedicated test added in Phase 2.
    UploadEndpoint(
        "/api/v1/reactions/enumerate_library_file",
        form={"reaction_smarts": "[C:1]Cl>>[C:1]O", "max_products": "5"},
        file_field="reactants_0",
        filename="reactants_0.smi",
        file_body="CCCl chloroethane\nCCCCl chloropropane\n",
        slow=True,
    ),
    # Phase 2 fixed parse_rxn_file to use AllChem.ReactionFromRxnBlock.
    UploadEndpoint(
        "/api/v1/reactions/upload_rxn",
        filename="reaction.rxn",
        file_body=RXN_FILE_BODY,
    ),
]


# --------------------------------------------------------------------------- #
# Static assets — references from base.html / app.js that must exist on disk
# --------------------------------------------------------------------------- #

STATIC_ASSETS: list[str] = [
    # Phase 3 design-token system. style.css was deleted in Phase 5.
    "static/css/tokens.css",
    "static/css/app.css",
    # Phase 4 JS stack — load order is significant (utils.js defines
    # window.NuGenUtils that the rest depend on).
    "static/js/utils.js",
    "static/js/formatters.js",
    "static/js/workspace.js",
    "static/js/results_card.js",
    "static/js/topbar_search.js",
    # Phase 5: legacy static/js/app.js was deleted (dead code; consumers
    # migrated to NuGenUtils + per-page modules). base_init.js is its
    # tiny replacement (lucide icon refresh + mobile menu wiring).
    "static/js/base_init.js",
    # Phase 5 per-page modules. Each pair (test_file_exists +
    # test_flask_serves) is auto-registered under tests/test_static_assets.py.
    "static/js/index.js",
    "static/js/descriptors.js",
    "static/js/fingerprints.js",
    "static/js/structure.js",
    "static/js/similarity.js",
    "static/js/reactions.js",
    "static/js/properties.js",
    "static/js/coordinates.js",
    "static/js/visualization.js",
    # Phase 7 polish: molecule picker modal on /structure.
    "static/js/structure_picker.js",
    "static/favicon.ico",
    "static/favicon.svg",
    "static/favicon-16x16.png",
    "static/favicon-32x32.png",
    "static/favicon-48x48.png",
    "static/favicon-64x64.png",
    "static/favicon-128x128.png",
    "static/favicon-192x192.png",
    "static/favicon-256x256.png",
    "static/favicon-512x512.png",
    "static/apple-touch-icon.png",
]


# --------------------------------------------------------------------------- #
# Convenience helpers
# --------------------------------------------------------------------------- #


def non_alias_endpoints() -> list[Endpoint]:
    """Endpoints that actually exercise a unique handler."""
    return [e for e in ALL_ENDPOINTS if not e.alias]


def fast_endpoints() -> list[Endpoint]:
    """Non-alias endpoints safe to run on every `pytest` invocation."""
    return [e for e in non_alias_endpoints() if not e.slow and not e.known_broken]


def fast_uploads() -> list[UploadEndpoint]:
    return [u for u in UPLOAD_ENDPOINTS if not u.slow and not u.known_broken]
