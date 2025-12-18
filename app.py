import streamlit as st
import pandas as pd
import numpy as np

from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdMolDescriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


DATA_PATH = "data/processed/lipophilicity_clean.csv"


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


def morgan_fp(smiles: str, radius: int = 2, nbits: int = 2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


@st.cache_resource
def load_data_and_index():
    df = pd.read_csv(DATA_PATH)
    df = df[df["smiles"].notna()].copy().reset_index(drop=True)

    fps = []
    keep = []
    for i, s in enumerate(df["smiles"].values):
        fp = morgan_fp(s)
        if fp is not None:
            fps.append(fp)
            keep.append(i)

    df = df.iloc[keep].reset_index(drop=True)
    return df, fps


@st.cache_resource
def train_model(df: pd.DataFrame):
    y = df["Y"].astype(float).values
    X = np.vstack([rdkit_features(s) for s in df["smiles"].values])

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    model = RandomForestRegressor(n_estimators=300, random_state=42, n_jobs=-1)
    model.fit(X_train, y_train)

    preds = model.predict(X_test)
    rmse = float(mean_squared_error(y_test, preds) ** 0.5)
    mae = float(mean_absolute_error(y_test, preds))
    r2 = float(r2_score(y_test, preds))
    return model, rmse, mae, r2


st.set_page_config(page_title="ChemiMLOps Demo", layout="centered")
st.title("🧪 ChemiMLOps — ADMET + Similarity Search")
st.caption("Lipophilicity prediction + top-K similarity search (RDKit + ML).")

df, fps = load_data_and_index()
model, rmse, mae, r2 = train_model(df)

with st.expander("Model metrics (baseline)", expanded=False):
    st.write(f"RMSE: **{rmse:.4f}** | MAE: **{mae:.4f}** | R²: **{r2:.4f}**")

smiles = st.text_input("SMILES", value="c1ccccc1")
k = st.slider("Top-K similar molecules", min_value=1, max_value=20, value=5)

col1, col2 = st.columns(2)

with col1:
    if st.button("Predict Lipophilicity"):
        X = rdkit_features(smiles)
        if X is None:
            st.error("Invalid SMILES. Example: CCO or c1ccccc1")
        else:
            pred = float(model.predict(X)[0])
            st.success(f"Predicted lipophilicity: **{pred:.4f}**")

with col2:
    if st.button("Find Similar Molecules"):
        qfp = morgan_fp(smiles)
        if qfp is None:
            st.error("Invalid SMILES for similarity search.")
        else:
            sims = DataStructs.BulkTanimotoSimilarity(qfp, fps)
            top_idx = sorted(range(len(sims)), key=lambda i: sims[i], reverse=True)[:k]
            results = []
            for i in top_idx:
                results.append(
                    {
                        "smiles": df.loc[i, "smiles"],
                        "inchikey": df.loc[i, "inchikey"],
                        "similarity": float(sims[i]),
                        "Y": float(df.loc[i, "Y"]),
                    }
                )
            st.subheader("Top Similar Molecules")
            st.dataframe(results, use_container_width=True)
