from rdkit import Chem
from rdkit.Chem import AllChem

class MoleculeReader:
    @staticmethod
    def read_from_smiles(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol)
        return ("mol", Chem.MolToMolBlock(mol), mol)

    @staticmethod
    def read_from_file(file_path):
        if file_path.endswith(".mol") or file_path.endswith(".sdf"):
            suppl = Chem.SDMolSupplier(file_path)
            mol = next((m for m in suppl if m is not None), None)
            if mol is None:
                raise ValueError("Failed to read molecule from file")
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.UFFOptimizeMolecule(mol)
            return ("mol", Chem.MolToMolBlock(mol), mol)
        elif file_path.endswith(".pdb"):
            with open(file_path, 'r') as f:
                return ("pdb", f.read(), None)
        else:
            raise ValueError("Unsupported file format")
