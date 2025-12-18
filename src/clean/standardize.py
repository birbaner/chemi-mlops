from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def canonicalize_smiles(smiles: str) -> str | None:
    """Return canonical SMILES or None if invalid."""
    if smiles is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    # sanitize + canonical SMILES
    return Chem.MolToSmiles(mol, canonical=True)


def smiles_to_inchikey(smiles: str) -> str | None:
    """Stable ID for deduplication (chemical data management)."""
    can = canonicalize_smiles(smiles)
    if can is None:
        return None
    mol = Chem.MolFromSmiles(can)
    return Chem.inchi.MolToInchiKey(mol)


def morgan_fp(smiles: str, radius: int = 2, nbits: int = 2048):
    """Fingerprint for similarity search."""
    can = canonicalize_smiles(smiles)
    if can is None:
        return None
    mol = Chem.MolFromSmiles(can)
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


if __name__ == "__main__":
    tests = ["CCO", "c1ccccc1", "not_a_smiles"]
    for s in tests:
        print(s, "=>", canonicalize_smiles(s), "=>", smiles_to_inchikey(s))
