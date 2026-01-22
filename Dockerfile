FROM dolfinx/dolfinx:stable

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

# install python deps (cache friendly)
COPY requirements.txt /app/requirements.txt
RUN python -m pip install --no-cache-dir -r /app/requirements.txt

# copy project files (optional for dev, useful for CI)
COPY . /app

# default to an interactive shell
CMD ["bash"]