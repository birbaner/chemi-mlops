FROM python:3.11-slim

# System dependencies for RDKit
RUN apt-get update && apt-get install -y \
    build-essential \
    libgl1 \
    libxrender1 \
    libxext6 \
    libsm6 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy project files
COPY . /app

# Install Python dependencies
RUN pip install --upgrade pip \
    && pip install -r requirements.txt

# Hugging Face Spaces expects port 7860
EXPOSE 7860

# IMPORTANT: run Streamlit on the correct port
CMD ["bash", "-lc", "streamlit run app.py --server.address 0.0.0.0 --server.port 7860 --server.headless true --browser.gatherUsageStats false"]
