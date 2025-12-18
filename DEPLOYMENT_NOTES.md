# DEPLOYMENT_NOTES.md

## Chemi Mlops Demo: End-to-End Deployment, Debugging, and Fixes

This document summarizes the process, root causes, and solutions for deploying the Chemi Mlops Demo (lipophilicity prediction + top-K similarity search) as a Docker-based Hugging Face Space. It is based on the detailed notes in the README and is intended as a standalone deployment reference.

---

### 1. Symptom
- Hugging Face Space status stuck at "Starting".
- Container logs in the web UI returned the HTML frontend, but container stdout was not visible.
- Space intermittently timed out during startup.

### 2. Root Causes
- Heavy data-loading and model training executed at module import time in `app.py`, blocking Streamlit from responding to readiness probes.
- Top-level import of the `tdc` package in `src/ingest/tdc_loader.py` (not listed in `requirements.txt` originally) caused inconsistent build/runtime behavior.
- Streamlit bind address and Docker HEALTHCHECK probing mismatches prevented the Hugging Face frontend from detecting readiness.

### 3. Fixes Applied

#### a. Lazy Import and Dependency Fixes
- Moved `from tdc.single_pred import ADME` inside the loader function in `src/ingest/tdc_loader.py` to avoid import-time failures.
- Added `tdc` to `requirements.txt`.

#### b. Avoid Long Work at Import Time
- Deferred heavy tasks (data loading and model training) to runtime inside a `st.spinner(...)` in `app.py` so Streamlit can start quickly and health probes can succeed.

#### c. Docker / Streamlit Runtime Adjustments
- Updated `Dockerfile` to bind Streamlit to `0.0.0.0` and use the `PORT` env var (fallback to `8501`).
- Removed or avoided unstable custom HEALTHCHECKs that probed the wrong host/port.

#### d. Temporary Debug Logging
- Added short startup logs and a `READY` marker to container logs for visibility during debugging. (Remove these when finished.)

### 4. Files Changed (High Level)
- `app.py` — moved data loading and training into a runtime spinner.
- `src/ingest/tdc_loader.py` — made `tdc` import lazy.
- `requirements.txt` — added `tdc`.
- `Dockerfile` — updated bind address, port usage, and removed debug HEALTHCHECK.

### 5. Validation Steps
- Fetched container run logs using the Hugging Face iframe logs endpoint and browser DevTools.
- Observed Streamlit startup messages and verified the app UI loaded and accepted input.
- Confirmed the app no longer performed heavy startup work during import and that the Hugging Face frontend detected readiness.

### 6. Best Practices & Next Steps
- Never perform expensive CPU/IO operations at import time in containerized apps.
- Keep `requirements.txt` complete and optionally pin versions.
- For explicit readiness probes, add a minimal HTTP endpoint that returns 200 only when ready; point Docker HEALTHCHECK to it.
- Remove development debug logging and HEALTHCHECK workarounds once stable.

### 7. Commit & Push Instructions

```powershell
git add README.md Dockerfile src/ingest/tdc_loader.py requirements.txt app.py
git commit -m "docs: add end-to-end deployment notes and changelog"
git push hf main
```

---

This document is derived from the deployment and debugging notes in the README. For further details, see the main README file.
