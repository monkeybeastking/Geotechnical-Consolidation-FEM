# --- Base image with Dolfinx preinstalled ---
FROM dolfinx/dolfinx:stable

# --- Optional environment variable for clarity ---
ARG PYTHON_ENV=my_env
ENV PYTHON_ENV=$PYTHON_ENV

# --- Working directory ---
WORKDIR /workspaces/app

# --- Copy your files ---
COPY requirements.txt ./requirements.txt

# --- Install extra Python tools you want ---
RUN pip install --no-cache-dir -r requirements.txt
