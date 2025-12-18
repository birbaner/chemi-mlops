from pathlib import Path
import pandas as pd
import numpy as np
import joblib

from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor


def rdkit_features(smiles: str) -> list[float]:
    mol = Chem.MolFromSmiles(smiles)
    # simple, strong baseline descriptors
    return [
        Descriptors.MolWt(mol),
        Descriptors.MolLogP(mol),
        Descriptors.TPSA(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        Descriptors.NumRotatableBonds(mol),
        Descriptors.RingCount(mol),
    ]


if __name__ == "__main__":
    df = pd.read_csv("data/processed/lipophilicity_clean.csv")

    # Target column in TDC is usually "Y"
    y = df["Y"].astype(float).values
    X = np.array([rdkit_features(s) for s in df["smiles"].values], dtype=float)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    model = RandomForestRegressor(
        n_estimators=300,
        random_state=42,
        n_jobs=-1
    )
    model.fit(X_train, y_train)

    preds = model.predict(X_test)

    rmse = mean_squared_error(y_test, preds) ** 0.5

    mae = mean_absolute_error(y_test, preds)
    r2 = r2_score(y_test, preds)

    print(f"RMSE: {rmse:.4f}")
    print(f"MAE : {mae:.4f}")
    print(f"R^2 : {r2:.4f}")

    Path("models").mkdir(exist_ok=True)
    joblib.dump(model, "models/lipophilicity_rf.joblib")
    print("Saved model to models/lipophilicity_rf.joblib")
