FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

COPY . /app

RUN pip install --upgrade pip && pip install -r requirements.txt

EXPOSE 7860

# Temporarily remove HEALTHCHECK to avoid premature container restarts
# HEALTHCHECK CMD curl --fail http://localhost:${PORT:-7860}/_stcore/health || exit 1

CMD ["sh", "-c", "PORT=${PORT:-${STREAMLIT_SERVER_PORT:-7860}}; echo \"ENV: PORT=$PORT STREAMLIT_SERVER_PORT=$STREAMLIT_SERVER_PORT\"; python - <<'PY'\nimport os,sys\nprint('PY ENV PORT=',repr(os.environ.get('PORT')) )\nprint('PY ENV STREAMLIT_SERVER_PORT=',repr(os.environ.get('STREAMLIT_SERVER_PORT')) )\nPY\n; echo \"Starting Streamlit on port $PORT\"; streamlit run app.py --server.address=0.0.0.0 --server.port=$PORT --server.headless=true --server.enableCORS=false --server.enableXsrfProtection=false --browser.gatherUsageStats=false"]
