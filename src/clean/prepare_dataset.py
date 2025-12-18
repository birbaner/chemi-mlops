from pathlib import Path
import pandas as pd
from src.clean.standardize import canonicalize_smiles, smiles_to_inchikey


def prepare_csv(in_csv: str, out_csv: str) -> None:
    df = pd.read_csv(in_csv)

    # TDC commonly uses: Drug (SMILES) and Y (label)
    if "Drug" not in df.columns:
        raise ValueError(f"Expected column 'Drug' in {in_csv}. Found: {df.columns.tolist()}")

    df["smiles_raw"] = df["Drug"].astype(str)
    df["smiles"] = df["smiles_raw"].apply(canonicalize_smiles)
    df["inchikey"] = df["smiles"].apply(lambda s: smiles_to_inchikey(s) if s else None)

    # Drop invalid SMILES
    df = df[df["smiles"].notna()].copy()

    # Deduplicate by inchikey (keeps first)
    df = df.drop_duplicates(subset=["inchikey"])

    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)

    print("Rows saved:", len(df))
    print("Saved to:", out_csv)


if __name__ == "__main__":
    raw = "data/raw/tdc_adme_Lipophilicity_AstraZeneca.csv"
    out = "data/processed/lipophilicity_clean.csv"
    prepare_csv(raw, out)
