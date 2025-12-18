---
title: Chemi Mlops Demo
emoji: 🧪
colorFrom: red
colorTo: blue
sdk: docker
app_port: 8501
app_file: app.py
tags:
  - streamlit
  - cheminformatics
  - drug-discovery
  - rdkit
  - machine-learning
  - interpretability
pinned: false
short_description: Lipophilicity prediction with SHAP and similarity
---

# 🧪 ChemiMLOps Demo

**Professional computational drug discovery platform** featuring lipophilicity prediction, molecular similarity search, and ML interpretability.

Built with RDKit, scikit-learn, SHAP, and Streamlit. Deployed as a Docker-based Hugging Face Space.

---

## 🚀 Features

### 🧪 Molecular Descriptor Calculation
Automated calculation of key molecular descriptors using RDKit:
- **MolWt**: Molecular weight
- **LogP**: Octanol-water partition coefficient (lipophilicity)
- **TPSA**: Topological polar surface area
- **NumHDonors/NumHAcceptors**: Hydrogen bond donors/acceptors
- **NumRotatableBonds**: Rotatable bonds count
- **RingCount**: Number of rings

**Usage:**
```bash
# Integrated into data pipeline
python src/clean/prepare_dataset.py

# Standalone CLI
python scripts/add_descriptors_cli.py input.csv output.csv
```

### 🤖 Model Interpretability & Explainability
Professional-grade ML interpretability using **SHAP** (SHapley Additive exPlanations):
- **Global Feature Importance**: Bar chart showing which molecular properties drive predictions across the entire dataset
- **Local Explanations**: SHAP values for individual predictions showing feature contributions
- **Uncertainty Quantification**: 95% prediction intervals using Random Forest ensemble variance
- **Interactive Visualization**: SHAP waterfall plots with color-coded positive/negative contributions

**Key capabilities:**
- Understand *why* the model makes specific predictions
- Identify which molecular features increase/decrease lipophilicity
- Quantify prediction confidence with uncertainty estimates
- Communicate model decisions to chemists and stakeholders

### 🔍 Molecular Similarity Search
Fast similarity search using Morgan fingerprints and Tanimoto similarity:
- Top-K nearest neighbors for any query SMILES
- Pre-computed fingerprint index for 4,200+ compounds
- Returns similarity scores and lipophilicity values

---

## 📁 Project Structure

```
chemi-mlops/
├── app.py                          # Main Streamlit application
├── Dockerfile                      # Docker configuration for HF Space
├── requirements.txt                # Python dependencies
├── data/
│   ├── processed/
│   │   └── lipophilicity_clean.csv # Cleaned dataset with descriptors
│   └── raw/                        # Raw data (not tracked)
├── src/
│   ├── clean/
│   │   ├── prepare_dataset.py      # Data preparation pipeline
│   │   ├── add_descriptors.py      # Molecular descriptor calculation
│   │   └── standardize.py          # SMILES canonicalization
│   ├── interpretability/
│   │   └── shap_explainer.py       # SHAP-based model interpretation
│   ├── serving/
│   │   ├── api.py                  # FastAPI endpoints
│   │   └── similarity.py           # Similarity search engine
│   └── train/
│       └── train_baseline.py       # Model training script
└── scripts/
    └── add_descriptors_cli.py      # CLI for descriptor enrichment
```

---

## 🛠️ Technical Stack

- **Cheminformatics**: RDKit for molecular descriptor calculation and fingerprinting
- **Machine Learning**: scikit-learn RandomForestRegressor (300 estimators)
- **Interpretability**: SHAP TreeExplainer for feature attribution
- **Visualization**: Matplotlib, Altair, Streamlit
- **Data**: Therapeutics Data Commons (TDC) Lipophilicity dataset (AstraZeneca)
- **Deployment**: Docker + Hugging Face Spaces

---

## 🎯 Use Cases

This project demonstrates professional capabilities for:

1. **Computational Drug Discovery**: Predict ADMET properties (lipophilicity) from molecular structure
2. **Cheminformatics Engineering**: Handle molecular data, compute descriptors, perform similarity searches
3. **Explainable AI in Chemistry**: Interpret black-box models for scientific decision-making
4. **MLOps for Life Sciences**: Deploy reproducible ML pipelines in containerized environments

---

## 📊 Model Performance

- **RMSE**: ~0.55 (lipophilicity units)
- **MAE**: ~0.42
- **R²**: ~0.75
- **Dataset**: 4,200 molecules from AstraZeneca (via TDC)
- **Features**: 7 RDKit molecular descriptors

---

## 🚢 Deployment

This section captures the investigation, root causes, and final fixes that resolved a persistent "Starting" state in the HF Space UI.

- **Symptom:** HF Space status stayed at "Starting" while the container logs shown in the web UI returned the HTML frontend; container stdout was not visible in the Logs panel and the Space intermittently timed out during startup.

- **Root causes discovered:**
  - Heavy data-loading and model training executed at module import time in `app.py`, blocking the Streamlit process from responding to readiness probes.
  - A top-level import of the `tdc` package in `src/ingest/tdc_loader.py` (and it not being listed in `requirements.txt` originally) caused inconsistent build/runtime behavior.
  - Streamlit bind address and Docker HEALTHCHECK probing mismatches caused the Hugging Face frontend to fail to detect the app as ready in some experiments.

### Fixes applied

1. Lazy import and dependency fixes
	- Moved `from tdc.single_pred import ADME` inside the loader function to avoid import-time failures when `tdc` is not yet installed or available.
	- Added `tdc` to `requirements.txt` so the package is installed during the build.

2. Avoid long work at import time
	- Deferred heavy tasks (data loading and baseline model training) to runtime inside a `st.spinner(...)` so Streamlit can start quickly and health probes can succeed.

3. Docker / Streamlit runtime adjustments
	- Updated `Dockerfile` to bind Streamlit to `0.0.0.0` and use the `PORT` env var (fallback to `8501`) so the platform can reach the app.
	- Removed or avoided unstable custom HEALTHCHECKs used during debugging that probed the wrong host/port.

4. Temporary debug logging
	- Added short startup logs and a `READY` marker while iterating on fixes to make readiness visible in the container logs. Trim/clean these when finished.

### Files changed (high level)

- `app.py` — moved data loading and training into a runtime spinner instead of running on import.
- `src/ingest/tdc_loader.py` — moved `tdc` import into the loader function (lazy import).
- `requirements.txt` — added `tdc`.
- `Dockerfile` — changed `--server.address` to `0.0.0.0`, used `PORT` fallback `8501`, removed debug HEALTHCHECK.

### How validation was performed

- Fetched container run logs using the HF iframe logs endpoint (`?logs=container`) and via browser DevTools Network → Response when the endpoint returned HTML.
- Observed Streamlit startup messages and verified the app UI (iframe) loaded and accepted input (SMILES → similarity search + model metrics).
- Confirmed the app no longer performed heavy startup work during import and that the HF frontend detected the app as running when bound to `0.0.0.0`.

## 🚢 Deployment

### Hugging Face Space

Live demo: `https://huggingface.co/spaces/rb757/chemi-mlops-demo`

The app is deployed as a Docker container on Hugging Face Spaces. Every push to the `hf` remote triggers an automatic rebuild.

### Local Development

```bash
# Clone repository
git clone https://github.com/birbaner/chemi-mlops.git
cd chemi-mlops

# Install dependencies
pip install -r requirements.txt

# Run Streamlit app
streamlit run app.py
```

### Deployment Notes

For complete deployment debugging and troubleshooting, see [DEPLOYMENT_NOTES.md](DEPLOYMENT_NOTES.md).

**Key deployment fixes:**
- Lazy imports to avoid startup failures
- Deferred heavy workloads to runtime (inside `st.spinner`)
- Proper Docker networking (bind to `0.0.0.0`)
- Clean YAML frontmatter in README for HF Space configuration

---

## 📚 Best Practices Demonstrated

- ✅ **Modular Code Structure**: Separate modules for data, training, serving, and interpretability
- ✅ **Reproducible Pipelines**: Automated data preparation and descriptor calculation
- ✅ **Explainable AI**: SHAP integration for transparent ML predictions
- ✅ **Uncertainty Quantification**: Prediction intervals from ensemble methods
- ✅ **Production Deployment**: Dockerized app with proper networking and health checks
- ✅ **Version Control**: Git workflow with separate remotes for GitHub and Hugging Face

---

## 🔄 Git Workflow

```bash
# Stage and commit changes
git add .
git commit -m "feat: your feature description"

# Push to GitHub
git push origin main

# Push to Hugging Face Space (triggers rebuild)
git push hf main
```

---

## 📈 Future Extensions

- [ ] Multi-property ADMET prediction (solubility, permeability, toxicity)
- [ ] 2D molecular structure visualization
- [ ] Substructure and scaffold search
- [ ] Chemical space exploration (t-SNE/UMAP)
- [ ] Model versioning and experiment tracking (MLflow)
- [ ] API endpoints for programmatic access

---

## 📄 License

MIT License - see LICENSE file for details.

---

## 🙏 Acknowledgments

- **Therapeutics Data Commons (TDC)** for the Lipophilicity dataset
- **RDKit** for cheminformatics toolkit
- **SHAP** for model interpretability
- **Hugging Face** for Spaces deployment platform
