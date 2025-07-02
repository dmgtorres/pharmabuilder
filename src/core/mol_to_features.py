from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import os

FEATURE_TYPE_MAP = {
    "Donor": "HBD",
    "SingleAtomDonor": "HBD",
    "Acceptor": "HBA",
    "SingleAtomAcceptor": "HBA",
    "Aromatic": "AR",
    "Hydrophobe": "H",
    "LumpedHydrophobe": "H",
    "PosIonizable": "PI",
    "NegIonizable": "NI",
    "ZnBinder": "ZnB",
    "Imidazole": "AR",
    "Arom5": "AR"
}

def get_structured_features(mol):
    from rdkit.Chem import rdMolTransforms

    fdef = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
    factory = ChemicalFeatures.BuildFeatureFactory(fdef)

    feats = factory.GetFeaturesForMol(mol)
    structured = []
    ring_info = mol.GetRingInfo()

    for feat in feats:
        ftype = feat.GetType()
        tag = FEATURE_TYPE_MAP.get(ftype, "Unknown")
        atom_ids = list(feat.GetAtomIds())

        if tag == "Unknown":
            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in atom_ids):
                if any(set(atom_ids) == set(ring) for ring in ring_info.AtomRings()):
                    tag = "AR"
                else:
                    continue
            else:
                continue

        if tag == "H":
            valid_atoms = []
            for i in atom_ids:
                atom = mol.GetAtomWithIdx(i)
                if atom.GetAtomicNum() != 6:
                    valid_atoms.append(i)
                    continue
                if atom.GetDegree() >= 3:
                    if any(nbr.GetAtomicNum() in (7, 8) for nbr in atom.GetNeighbors()):
                        continue
                valid_atoms.append(i)

            if not valid_atoms:
                continue
            atom_ids = valid_atoms

        conf = mol.GetConformer()
        points = [conf.GetAtomPosition(i) for i in atom_ids]
        if not points:
            continue
        x = sum(p.x for p in points) / len(points)
        y = sum(p.y for p in points) / len(points)
        z = sum(p.z for p in points) / len(points)

        structured.append({
            "functional_group": feat.GetFamily(),
            "pharmacophore_features": [tag],
            "atom_ids": atom_ids,
            "centroid": (round(x, 3), round(y, 3), round(z, 3))
        })

    # Aromatic nitrogen fallback (manual HBA)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        if not atom.GetIsAromatic():
            continue
        if atom.GetTotalNumHs() > 0:
            continue
        if atom.GetHybridization().name != "SP2":
            continue
        if any(atom.GetIdx() in f["atom_ids"] for f in structured if "HBA" in f["pharmacophore_features"]):
            continue
        conf = mol.GetConformer()
        pos = conf.GetAtomPosition(atom.GetIdx())
        structured.append({
            "functional_group": "AromaticNitrogen",
            "pharmacophore_features": ["HBA"],
            "atom_ids": [atom.GetIdx()],
            "centroid": (round(pos.x, 3), round(pos.y, 3), round(pos.z, 3))
        })

    # Aromatic oxygen fallback (manual HBA)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        if not atom.GetIsAromatic():
            continue
        if atom.GetTotalNumHs() > 0:
            continue
        if atom.GetHybridization().name != "SP2":
            continue
        if any(atom.GetIdx() in f["atom_ids"] for f in structured if "HBA" in f["pharmacophore_features"]):
            continue
        conf = mol.GetConformer()
        pos = conf.GetAtomPosition(atom.GetIdx())
        structured.append({
            "functional_group": "AromaticOxygen",
            "pharmacophore_features": ["HBA"],
            "atom_ids": [atom.GetIdx()],
            "centroid": (round(pos.x, 3), round(pos.y, 3), round(pos.z, 3))
        })

    # Methyl (CH3) fallback as hydrophobe
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        h_count = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1)
        if h_count == 3:
            if any(atom.GetIdx() in f["atom_ids"] for f in structured if "H" in f["pharmacophore_features"]):
                continue
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            structured.append({
                "functional_group": "Hydrophobic",
                "pharmacophore_features": ["H"],
                "atom_ids": [atom.GetIdx()],
                "centroid": (round(pos.x, 3), round(pos.y, 3), round(pos.z, 3))
            })

    return structured
