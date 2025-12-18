from pathlib import Path
import pandas as pd
from tdc.single_pred import ADME


def load_and_save_adme_dataset(name: str, out_dir: str = "data/raw") -> Path:
    """
    Downloads an ADME dataset from TDC and saves it locally as CSV.

    Typical columns include:
      - Drug (SMILES)
      - Y (label/target)
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    data = ADME(name=name)
    df = data.get_data()

    csv_path = out_path / f"tdc_adme_{name}.csv"
    df.to_csv(csv_path, index=False)
    return csv_path


if __name__ == "__main__":
    # A common starter dataset in ADME (you can change later)
    dataset_name = "Lipophilicity_AstraZeneca"
    saved = load_and_save_adme_dataset(dataset_name)
    print(f"Saved dataset to: {saved}")
