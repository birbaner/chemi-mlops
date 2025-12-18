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

CMD ["sh", "-c", "streamlit run app.py --server.address=0.0.0.0 --server.port=${PORT:-7860} --server.headless=true --server.enableCORS=false --server.enableXsrfProtection=false --browser.gatherUsageStats=false"]
