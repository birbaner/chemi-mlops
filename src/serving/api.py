from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field
import joblib
import numpy as np

from rdkit import Chem
from rdkit.Chem import Descriptors

from src.serving.similarity import SimilarityIndex

app = FastAPI(title="ChemiMLOps API", version="0.1.0")

MODEL_PATH = "models/lipophilicity_rf.joblib"
DATA_PATH = "data/processed/lipophilicity_clean.csv"

model = joblib.load(MODEL_PATH)
SIM_INDEX: SimilarityIndex | None = None


class PredictRequest(BaseModel):
    smiles: str = Field(..., example="CCO")  # ethanol example


def rdkit_features(smiles: str) -> np.ndarray | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    feats = [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.TPSA(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.RingCount(mol),
    ]
    return np.array(feats, dtype=float).reshape(1, -1)


@app.on_event("startup")
def load_similarity_index():
    """Load similarity search index once when the API starts."""
    global SIM_INDEX
    SIM_INDEX = SimilarityIndex(DATA_PATH)


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/predict/lipophilicity")
def predict(req: PredictRequest):
    X = rdkit_features(req.smiles)
    if X is None:
        raise HTTPException(
            status_code=400,
            detail="Invalid SMILES. Example: 'CCO' or 'c1ccccc1'",
        )

    pred = float(model.predict(X)[0])
    return {"smiles": req.smiles, "lipophilicity_pred": pred}


@app.get("/similarity/topk")
def similarity_topk(smiles: str, k: int = 5):
    if SIM_INDEX is None:
        raise HTTPException(status_code=503, detail="Similarity index not loaded yet")

    results = SIM_INDEX.topk(smiles, k=k)
    if results is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES for similarity search")

    return {"query": smiles, "k": k, "results": results}
