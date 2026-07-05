# PharmInform container image.
#
# Uses python:3.12-slim rather than matching this project's development
# Python (3.14) -- 3.12 has mature, widely-available wheels/apt packages
# for every native dependency below (cyvcf2, weasyprint's Pango/cairo
# stack), which was not a given for 3.14 at the time this was written.
FROM python:3.12-slim

# Native libraries required by weasyprint (PDF export) and cyvcf2 (htslib
# build deps). Confirmed necessary from this project's own development
# experience: weasyprint failed to import with a bare `pip install` until
# the Pango/GObject stack was present (see report.py's docstring for the
# equivalent macOS/Homebrew issue).
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libpango-1.0-0 \
    libpangoft2-1.0-0 \
    libgdk-pixbuf-2.0-0 \
    libcairo2 \
    libffi-dev \
    shared-mime-info \
    fonts-liberation \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY pyproject.toml README.md ./
COPY src ./src
COPY templates ./templates
COPY streamlit_app ./streamlit_app

RUN pip install --no-cache-dir .

# Aldy is an optional extra (heavy dependency, requires an ILP solver) --
# not installed by default. Uncomment to include it:
# RUN pip install --no-cache-dir .[aldy]

EXPOSE 8000 8501

# Default: serve the FastAPI app. Override the command to run the CLI or
# the Streamlit demo instead, e.g.:
#   docker run pharminform pharminform report --vcf /data/sample.vcf --out /data/report.html
#   docker run -p 8501:8501 pharminform streamlit run streamlit_app/app.py --server.address 0.0.0.0
CMD ["uvicorn", "pharminform.api:app", "--host", "0.0.0.0", "--port", "8000"]
