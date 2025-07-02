from rdkit import Chem
from rdkit.Chem import AllChem

class MoleculeStore:
    def __init__(self):
        self._mol = None
        self._mol_block = None
        self._mol_type = None

    def load_from_smiles(self, smiles: str):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)
        self._mol = mol
        self._mol_block = Chem.MolToMolBlock(mol)
        self._mol_type = "mol"

    def load_from_file(self, file_path: str):
        if file_path.endswith(".mol") or file_path.endswith(".sdf"):
            mol = Chem.MolFromMolFile(file_path, sanitize=True)
        elif file_path.endswith(".pdb"):
            mol = Chem.MolFromPDBFile(file_path, sanitize=True)
        else:
            raise ValueError("Unsupported file format")

        if mol is None:
            raise ValueError("Could not parse file")

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)
        self._mol = mol
        self._mol_block = Chem.MolToMolBlock(mol)
        self._mol_type = "mol" if file_path.endswith(".mol") or file_path.endswith(".sdf") else "pdb"

    def get_rdkit_mol(self):
        return self._mol

    def get_viewable_data(self):
        if self._mol and self._mol_block and self._mol_type:
            return self._mol_type, self._mol_block, self._mol
        raise ValueError("Molecule not fully initialized.")
