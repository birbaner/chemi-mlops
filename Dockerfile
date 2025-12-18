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

HEALTHCHECK --start-period=60s --interval=10s --timeout=5s --retries=5 CMD curl --fail http://localhost:7860/_stcore/health || exit 1

CMD ["sh", "-c", "PORT=${PORT:-7860}; sleep 20 && streamlit run app.py --server.address=0.0.0.0 --server.port=$PORT --server.headless=true --server.enableCORS=false --server.enableXsrfProtection=false --browser.gatherUsageStats=false"]
