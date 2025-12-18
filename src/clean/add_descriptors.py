import sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


def compute_descriptors(smiles: str) -> dict:
    """Compute a set of molecular descriptors for a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    desc = {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "TPSA": Descriptors.TPSA(mol),
    }
    return desc


def add_descriptors_to_csv(in_csv: str, out_csv: str) -> None:
    df = pd.read_csv(in_csv)
    descs = df["smiles"].apply(compute_descriptors)
    desc_df = pd.DataFrame(list(descs))
    df = pd.concat([df, desc_df], axis=1)
    df.to_csv(out_csv, index=False)
