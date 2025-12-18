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

# Bind to 0.0.0.0 so the HF frontend can reach Streamlit. Use PORT env var (fallback 8501).
CMD ["sh", "-c", "PORT=${PORT:-8501}; exec streamlit run app.py --server.address=0.0.0.0 --server.port=$PORT --server.headless=true --server.enableCORS=false --server.enableXsrfProtection=false --browser.gatherUsageStats=false"]
