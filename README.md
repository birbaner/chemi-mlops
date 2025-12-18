# ChemiMLOps — ADMET (Lipophilicity) Prediction + Similarity Search

End-to-end cheminformatics ML project:
- Ingests ADMET data (TDC)
- Cleans + canonicalizes SMILES (RDKit) and generates InChIKey
- Trains an ML baseline model (RandomForest)
- Serves predictions via FastAPI
- Provides similarity search (Morgan fingerprints, Tanimoto)
- Includes Streamlit web demo UI

## Tech Stack
Python, RDKit, PyTDC, scikit-learn, FastAPI, Uvicorn, Streamlit

## Run locally
```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt

python -m src.ingest.tdc_loader
python -m src.clean.prepare_dataset
python -m src.train.train_baseline

uvicorn src.serving.api:app --reload
streamlit run app/streamlit_app.py
