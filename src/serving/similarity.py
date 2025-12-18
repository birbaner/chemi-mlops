import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors


def morgan_fp(smiles: str, radius: int = 2, nbits: int = 2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


class SimilarityIndex:
    def __init__(self, csv_path: str):
        self.df = pd.read_csv(csv_path)
        self.df = self.df[self.df["smiles"].notna()].copy()

        fps = []
        keep_rows = []
        for i, s in enumerate(self.df["smiles"].values):
            fp = morgan_fp(s)
            if fp is not None:
                fps.append(fp)
                keep_rows.append(i)

        self.df = self.df.iloc[keep_rows].reset_index(drop=True)
        self.fps = fps

    def topk(self, query_smiles: str, k: int = 5):
        qfp = morgan_fp(query_smiles)
        if qfp is None:
            return None

        sims = DataStructs.BulkTanimotoSimilarity(qfp, self.fps)
        top_idx = sorted(range(len(sims)), key=lambda i: sims[i], reverse=True)[:k]

        out = []
        for i in top_idx:
            out.append({
                "smiles": self.df.loc[i, "smiles"],
                "inchikey": self.df.loc[i, "inchikey"],
                "similarity": float(sims[i]),
                "Y": float(self.df.loc[i, "Y"]),
            })
        return out
