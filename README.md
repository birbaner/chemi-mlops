title: Chemi Mlops Demo
emoji: 🚀
colorFrom: red
colorTo: red
sdk: docker
app_port: 8501
tags:
pinned: false
short_description: Streamlit template space

title: Chemi Mlops Demo
emoji: 🚀
colorFrom: red
colorTo: red
sdk: docker
app_port: 8501
tags:
pinned: false
short_description: Streamlit template space

# Chemi Mlops Demo

---
title: Chemi Mlops Demo
emoji: 🚀
colorFrom: red
colorTo: red
sdk: docker
app_port: 8501
tags: ["streamlit"]
pinned: false
short_description: Streamlit template space
---

# Chemi Mlops Demo

Lipophilicity prediction + top-K similarity search (RDKit + ML).

This repository contains a Streamlit app and supporting code used to build a Docker-based Hugging Face Space. The sections below document the end-to-end debugging and fixes applied during deployment, and provide concise next steps and push instructions.

## End-to-End: Space deployment, debugging and fixes

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

## Changelog (concise)

- Initial: Streamlit template + baseline model training logic.
- Added `tdc` to `requirements.txt`; made import lazy to avoid import-time failures.
- Deferred heavy data processing/training into runtime spinner in `app.py`.
- Updated `Dockerfile` to make Streamlit reachable by HF (bind address and PORT usage).
- Temporary debug logging used while iterating; cleaned up in final commit.

## Cheminformatics Extensions

### Molecular Descriptor Calculation
The project now includes automated molecular descriptor calculation using RDKit. The following descriptors are computed for each molecule:
- **MolWt**: Molecular weight
- **LogP**: Octanol-water partition coefficient
- **NumHDonors**: Number of hydrogen bond donors
- **NumHAcceptors**: Number of hydrogen bond acceptors
- **TPSA**: Topological polar surface area

**Usage:**
1. Integrated into data pipeline: Run `python src/clean/prepare_dataset.py` to automatically enrich your dataset with descriptors.
2. Standalone CLI: `python scripts/add_descriptors_cli.py input.csv output.csv`

This extension demonstrates professional-grade cheminformatics capabilities suitable for computational drug discovery workflows.

## Best practices & next steps

- Never perform expensive CPU/IO operations at import time in apps that run under container orchestration or platform health probes.
- Keep `requirements.txt` complete (and optionally pin versions) for reproducible builds.
- If you need an explicit readiness probe, add a minimal HTTP endpoint (Flask/FastAPI) that returns 200 only when all required resources are ready; point Docker HEALTHCHECK to it.
- Remove development debug logging and HEALTHCHECK workarounds once the platform consistently detects readiness.

## Commit & push instructions

Run these commands locally to commit and push the README and other final changes to the `hf` remote:

```powershell
git add README.md Dockerfile src/ingest/tdc_loader.py requirements.txt app.py
git commit -m "docs: add end-to-end deployment notes and changelog"
git push hf main
```

If you prefer, I can craft a smaller commit message or split the commit into implementation vs docs steps — tell me which you prefer.

---
