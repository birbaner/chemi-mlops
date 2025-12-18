# DEPLOYMENT_NOTES.md

## ChemiMLOps Demo: End-to-End Deployment, Debugging, and Fixes

This document summarizes the complete deployment process, root causes, and solutions for deploying the ChemiMLOps Demo as a Docker-based Hugging Face Space.

---

## 📦 Project Overview

**Live Demo**: `https://huggingface.co/spaces/rb757/chemi-mlops-demo`  
**GitHub**: `https://github.com/birbaner/chemi-mlops`

**Tech Stack**:
- RDKit for cheminformatics
- scikit-learn for ML
- SHAP for interpretability
- Streamlit for UI
- Docker for deployment
- Hugging Face Spaces for hosting

---

## 🚨 Initial Deployment Issues

### 1. Symptom
- Hugging Face Space status stuck at "Starting"
- Container logs in the web UI returned the HTML frontend instead of stdout
- Space intermittently timed out during startup
- App failed to respond to health probes

### 2. Root Causes

#### a. Import-Time Heavy Operations
- Data loading and model training executed at module import time in `app.py`
- Blocked Streamlit from starting and responding to readiness probes
- HF platform couldn't detect when app was ready

#### b. Missing Dependencies
- Top-level import of `tdc` package in `src/ingest/tdc_loader.py`
- Package not listed in `requirements.txt`
- Caused inconsistent build/runtime behavior

#### c. Network Configuration Issues
- Streamlit bind address and Docker HEALTHCHECK mismatches
- App listening on localhost instead of `0.0.0.0`
- HF frontend couldn't reach the app

#### d. Invalid YAML Metadata
- Duplicate and malformed YAML frontmatter blocks in README
- `tags` field not properly formatted as array
- Caused repo card validation failures

---

## ✅ Solutions Implemented

### 1. Lazy Import and Dependency Fixes

**Problem**: Missing `tdc` dependency causing import failures

**Solution**:
```python
# src/ingest/tdc_loader.py
def load_lipophilicity():
    from tdc.single_pred import ADME  # Lazy import
    data = ADME(name="Lipophilicity_AstraZeneca")
    return data.get_data()
```

**Action**: Added `tdc` to `requirements.txt`

### 2. Deferred Heavy Workloads

**Problem**: Import-time model training blocking Streamlit startup

**Solution**:
```python
# app.py
with st.spinner("Loading data and training model..."):
    df, fps = load_data_and_index()
    model, interpreter, rmse, mae, r2 = train_model(df)
```

**Result**: Streamlit can start quickly and respond to health probes

### 3. Docker Networking

**Problem**: App not reachable from HF platform

**Solution**:
```dockerfile
# Dockerfile
CMD ["sh", "-c", "PORT=${PORT:-8501}; exec streamlit run app.py --server.address=0.0.0.0 --server.port=$PORT ...]
```

**Key changes**:
- Bind to `0.0.0.0` instead of `127.0.0.1`
- Use `PORT` environment variable with fallback to `8501`
- Enable proper CORS and XSRF settings for iframe embedding

### 4. Clean YAML Frontmatter

**Problem**: Repo card validation failures

**Solution**:
```yaml
---
title: Chemi Mlops Demo
emoji: 🧪
sdk: docker
app_port: 8501
app_file: app.py
tags:
  - streamlit
  - cheminformatics
---
```

**Result**: Single, properly formatted YAML block at top of README

---

## 🔧 Module Import Path Fixes

**Problem**: `ModuleNotFoundError: No module named 'src'` when running scripts directly

**Solution**: Added project root to `sys.path` in all entry point scripts:
```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
```

**Files updated**:
- `src/clean/prepare_dataset.py`
- `src/clean/add_descriptors.py`
- `scripts/add_descriptors_cli.py`

---

## 🎨 Feature Additions

### Phase 1: Molecular Descriptors (Cheminformatics)
- Created `src/clean/add_descriptors.py` for RDKit descriptor calculation
- Integrated into data preparation pipeline
- Added standalone CLI: `scripts/add_descriptors_cli.py`

**Descriptors computed**:
- MolWt, LogP, TPSA, NumHDonors, NumHAcceptors, NumRotatableBonds, RingCount

### Phase 2: Model Interpretability (ML/AI)
- Created `src/interpretability/shap_explainer.py` with `ModelInterpreter` class
- Integrated SHAP TreeExplainer for feature attribution
- Added uncertainty quantification using RF ensemble predictions
- Built interactive SHAP waterfall plots

**Capabilities**:
- Global feature importance visualization
- Local SHAP explanations for individual predictions
- 95% prediction intervals
- Color-coded contribution tables

---

## 🐛 Bug Fixes

### TypeError in SHAP Waterfall Plot
**Error**: `TypeError: unsupported format string passed to numpy.ndarray.__format__`

**Cause**: `base_value` from SHAP explainer was numpy array, not scalar

**Fix**:
```python
base_value = float(self.explainer.expected_value)
```

---

## 📊 Files Modified (Summary)

### Core Application
- `app.py` — Added SHAP integration, uncertainty quantification, enhanced UI
- `requirements.txt` — Added `shap`, `matplotlib`

### Data Pipeline
- `src/clean/prepare_dataset.py` — Integrated descriptor calculation
- `src/clean/add_descriptors.py` — New module for molecular descriptors
- `scripts/add_descriptors_cli.py` — CLI wrapper

### Interpretability
- `src/interpretability/__init__.py` — Package init
- `src/interpretability/shap_explainer.py` — SHAP-based model interpretation

### Configuration
- `README.md` — Complete rewrite with professional structure
- `DEPLOYMENT_NOTES.md` — This document
- `.gitignore` — Enhanced security patterns, allow processed data
- `Dockerfile` — Network and port configuration

---

## ✅ Validation & Testing

### Manual Testing
1. Run locally: `streamlit run app.py`
2. Test prediction with SMILES input (e.g., `CCO`, `c1ccccc1`)
3. Verify SHAP explanations appear
4. Check prediction intervals display correctly
5. Test similarity search functionality

### Deployment Validation
1. Push to HF: `git push hf main`
2. Monitor build logs in HF Space UI
3. Wait for "Running" status
4. Test all features in deployed app
5. Check for console errors in browser DevTools

---

## 📈 Best Practices Applied

1. **Never perform expensive I/O at import time** — Defer to runtime with spinners
2. **Complete dependency management** — All packages in `requirements.txt`
3. **Proper networking for containers** — Bind to `0.0.0.0`, use env vars
4. **Clean configuration** — Valid YAML, no duplicates
5. **Modular code structure** — Separate concerns (data, train, serve, interpret)
6. **Error handling** — Graceful degradation, informative error messages
7. **Security** — Comprehensive `.gitignore`, no secrets in repo
8. **Documentation** — README for users, DEPLOYMENT_NOTES for developers

---

## 🔄 Git Workflow

```bash
# Stage changes
git add <files>

# Commit with conventional message
git commit -m "feat: add feature X"

# Push to GitHub
git push origin main

# Push to Hugging Face (triggers rebuild)
git push hf main
```

---

## 🚀 Future Considerations

- Pin dependency versions for reproducibility
- Add explicit health check endpoint (`/health`)
- Implement logging with structured output
- Add unit tests for core modules
- Create CI/CD pipeline for automated testing
- Monitor performance metrics (inference time, memory usage)

---

**Last Updated**: December 18, 2025  
**Status**: ✅ Deployed and running successfully
