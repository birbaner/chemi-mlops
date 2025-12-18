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

CMD ["sh", "-c", "PORT=${PORT:-${STREAMLIT_SERVER_PORT:-7860}}; echo \"ENV: PORT=$PORT STREAMLIT_SERVER_PORT=$STREAMLIT_SERVER_PORT\"; python -c \"import os; print('PY ENV PORT=', repr(os.environ.get('PORT'))); print('PY ENV STREAMLIT_SERVER_PORT=', repr(os.environ.get('STREAMLIT_SERVER_PORT')))\"; echo \"Starting Streamlit on port $PORT\"; streamlit run app.py --server.address=0.0.0.0 --server.port=$PORT --server.headless=true --server.enableCORS=false --server.enableXsrfProtection=false --browser.gatherUsageStats=false & while ! curl -sSf http://127.0.0.1:$PORT/_stcore/health >/dev/null 2>&1; do echo \"waiting for /_stcore/health\"; sleep 1; done; echo \"health ok\"; wait"]
